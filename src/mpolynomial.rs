use crate::field::Field;
use crate::polynomial::Polynomial;
use crate::FieldElement;
use serde::Serialize;
use std::collections::HashMap;
use std::rc::Rc;

#[derive(Clone, Serialize)]
pub struct MPolynomial<T: FieldElement> {
    field: Rc<Field<T>>,
    exp_map: HashMap<Vec<u32>, T>,
}

impl<T: FieldElement> MPolynomial<T> {
    pub fn new(field: &Rc<Field<T>>) -> MPolynomial<T> {
        MPolynomial {
            field: Rc::clone(field),
            exp_map: HashMap::new(),
        }
    }

    pub fn from_map(map: &HashMap<Vec<u32>, T>, field: &Rc<Field<T>>) -> MPolynomial<T> {
        let mut m = HashMap::new();
        for (k, v) in map {
            m.insert(k.clone(), v.clone());
        }
        MPolynomial {
            field: Rc::clone(field),
            exp_map: m,
        }
    }

    pub fn exps(&self) -> &HashMap<Vec<u32>, T> {
        &self.exp_map
    }

    // only works if both polynomials are trimmed
    pub fn is_equal(&self, p: &MPolynomial<T>) -> bool {
        if self.exp_map.len() != p.exps().len() {
            return false;
        }
        for (k, v) in p.exps() {
            if let Some(val) = self.exp_map.get(k) {
                if val != v {
                    return false;
                }
            } else {
                return false;
            }
        }
        true
    }

    pub fn trim(&mut self) {
        let zero = T::zero();
        let mut to_remove = Vec::new();
        for (k, v) in self.exp_map.iter() {
            if v == &zero {
                to_remove.push(k.clone());
            }
        }
        for k in to_remove {
            self.exp_map.remove(&k);
        }
    }

    fn trim_vars(v: &Vec<u32>) -> Vec<u32> {
        for i in (0..v.len()).rev() {
            if v[i] != 0 {
                let mut out: Vec<u32> = vec![0; i + 1];
                out.copy_from_slice(&v[0..(i + 1)]);
                return out;
            }
        }
        vec![0]
    }

    pub fn term(&mut self, coef: &T, exps: &Vec<u32>) -> &Self {
        let zero = T::zero();
        let t = Self::trim_vars(exps);
        let existing = if let Some(v) = self.exp_map.get(&t) {
            v
        } else {
            &zero
        };
        self.exp_map.insert(t, self.field.add(existing, coef));
        self.trim();
        self
    }

    pub fn add(&mut self, poly: &MPolynomial<T>) -> &Self {
        for (k, v) in poly.exps().iter() {
            self.term(v, k);
        }
        self
    }

    pub fn sub(&mut self, poly: &MPolynomial<T>) -> &Self {
        for (k, v) in poly.exps().iter() {
            self.term(&self.field.neg(v), k);
        }
        self
    }

    pub fn mul_scalar(&mut self, val: &T) -> &Self {
        for (_exps, coef) in self.exp_map.iter_mut() {
            *coef = self.field.mul(coef, val);
        }
        self
    }

    pub fn mul(&mut self, poly: &MPolynomial<T>) -> &Self {
        let mut new_exp = HashMap::new();
        let zero = T::zero();
        for (exps1, coef1) in poly.exps() {
            for (exps2, coef2) in self.exps() {
                let mut final_exps: Vec<u32> = Vec::new();
                for i in 0..std::cmp::max(exps1.len(), exps2.len()) {
                    let e1 = exps1.get(i).unwrap_or(&0);
                    let e2 = exps2.get(i).unwrap_or(&0);
                    final_exps.push(e1 + e2);
                }
                let e = Self::trim_vars(&final_exps);
                let coef = self.field.add(
                    new_exp.get(&e).unwrap_or(&zero),
                    &self.field.mul(coef1, coef2),
                );
                new_exp.insert(e, coef);
            }
        }
        self.exp_map = new_exp;
        self
    }

    pub fn eval(&self, points: &[T]) -> T {
        let mut out = T::zero();
        for (exps, coef) in self.exps() {
            let mut inter = coef.clone();
            for i in 0..exps.len() {
                if exps[i] == 0 {
                    continue;
                }
                inter = self
                    .field
                    .lmul(&inter, &self.field.exp(&points[i], &T::from_u32(exps[i])));
            }
            out = self.field.add(&out, &inter);
        }
        out
    }

    pub fn eval_symbolic(&self, polys: &[Polynomial<T>]) -> Polynomial<T> {
        let mut out = Polynomial::new(&self.field);
        let mut degrees: Vec<u32> = Vec::new();
        for exps in self.exps().keys() {
            for (i, e) in exps.iter().enumerate() {
                if i >= degrees.len() {
                    degrees.resize(i + 1, 0);
                }
                if e + 1 > degrees[i] {
                    degrees[i] = e + 1;
                }
            }
        }
        let mut one = Polynomial::new(&self.field);
        one.term(&T::one(), 0);
        let mut power_map = Vec::new();
        for (i, d) in degrees.iter().enumerate() {
            let mut poly_powers = Vec::new();
            poly_powers.push(one.clone());
            for j in 1..usize::try_from(*d).unwrap() {
                poly_powers.push(Polynomial::mul_fft(
                    &poly_powers[j - 1],
                    &polys[i],
                    &self.field,
                ));
            }
            power_map.push(poly_powers);
        }
        for (exps, coef) in self.exps() {
            let mut inter = Polynomial::new(&self.field);
            inter.term(coef, 0);
            for i in 0..exps.len() {
                if exps[i] == 0 {
                    continue;
                }
                inter = Polynomial::mul_fft(
                    &inter,
                    &power_map[i][usize::try_from(exps[i]).unwrap()],
                    &self.field,
                );
            }
            out.add(&inter);
        }
        out
    }

    pub fn from_poly(poly: &Polynomial<T>) -> MPolynomial<T> {
        let mut out = MPolynomial::new(poly.field());
        for (exp, coef) in poly.coefs().iter().enumerate() {
            let exps: Vec<u32> = vec![u32::try_from(exp).unwrap()];
            out.term(coef, &exps);
        }
        out
    }

    pub fn variables(count: u32, field: &Rc<Field<T>>) -> Vec<MPolynomial<T>> {
        let mut out = Vec::new();
        for i in 0..count {
            let mut p = MPolynomial::new(field);
            let mut exps: Vec<u32> = vec![0; usize::try_from(i).unwrap()];
            exps.push(1);
            p.term(&T::one(), &exps);
            out.push(p);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;

    use crate::BigIntElement;

    use super::*;

    #[test]
    fn should_add_sub_multipolynomials() {
        let p = BigIntElement(BigInt::from(101));
        let g = BigIntElement(BigInt::from(0));
        let f = Rc::new(Field::new(p, g));

        // 4x + 2y^2 + 9
        let mut poly1 = MPolynomial::new(&f);
        poly1.term(&BigIntElement(BigInt::from(4)), &vec![1]);
        poly1.term(&BigIntElement(BigInt::from(2)), &vec![0, 2]);
        poly1.term(&BigIntElement(BigInt::from(9)), &vec![0]);

        // 2x^2 + y^2 + 99
        let mut poly2 = MPolynomial::new(&f);
        poly2.term(&BigIntElement(BigInt::from(2)), &vec![2]);
        poly2.term(&BigIntElement(BigInt::from(1)), &vec![0, 2]);
        poly2.term(&BigIntElement(BigInt::from(99)), &vec![0]);

        // 2x^2 + 4x + 3y^2 + 7
        let mut expected_add = MPolynomial::new(&f);
        expected_add.term(&BigIntElement(BigInt::from(2)), &vec![2]);
        expected_add.term(&BigIntElement(BigInt::from(4)), &vec![1]);
        expected_add.term(&BigIntElement(BigInt::from(3)), &vec![0, 2]);
        expected_add.term(&BigIntElement(BigInt::from(108)), &vec![0]);

        // -2x^2 + y^2 + 4x - 90
        let mut expected_sub = MPolynomial::new(&f);
        expected_sub.term(&BigIntElement(BigInt::from(4)), &vec![1]);
        expected_sub.term(&BigIntElement(BigInt::from(-2)), &vec![2]);
        expected_sub.term(&BigIntElement(BigInt::from(1)), &vec![0, 2]);
        expected_sub.term(&BigIntElement(BigInt::from(-90)), &vec![0]);

        assert!(poly1.clone().add(&poly2).is_equal(&expected_add));
        assert!(poly1.clone().sub(&poly2).is_equal(&expected_sub));
    }

    #[test]
    fn should_mul_multipolynomials() {
        let p = BigIntElement(BigInt::from(101));
        let g = BigIntElement(BigInt::from(0));
        let f = Rc::new(Field::new(p, g));

        // 4x + 2y^2 + 9
        let mut poly1 = MPolynomial::new(&f);
        poly1.term(&BigIntElement(BigInt::from(4)), &vec![1]);
        poly1.term(&BigIntElement(BigInt::from(2)), &vec![0, 2]);
        poly1.term(&BigIntElement(BigInt::from(9)), &vec![0]);

        // 2x^2 + y^2 + 99
        let mut poly2 = MPolynomial::new(&f);
        poly2.term(&BigIntElement(BigInt::from(2)), &vec![2]);
        poly2.term(&BigIntElement(BigInt::from(1)), &vec![0, 2]);
        poly2.term(&BigIntElement(BigInt::from(99)), &vec![0]);

        // 8x^3 + 4xy^2 + 396x + 4x^2y^2 + 2y^4 + 198y^2 + 18x^2 + 9y^2 + 891
        let mut expected = MPolynomial::new(&f);
        expected.term(&BigIntElement(BigInt::from(8)), &vec![3]);
        expected.term(&BigIntElement(BigInt::from(4)), &vec![1, 2]);
        expected.term(&BigIntElement(BigInt::from(396)), &vec![1]);
        expected.term(&BigIntElement(BigInt::from(4)), &vec![2, 2]);
        expected.term(&BigIntElement(BigInt::from(2)), &vec![0, 4]);
        expected.term(&BigIntElement(BigInt::from(198)), &vec![0, 2]);
        expected.term(&BigIntElement(BigInt::from(18)), &vec![2]);
        expected.term(&BigIntElement(BigInt::from(9)), &vec![0, 2]);
        expected.term(&BigIntElement(BigInt::from(891)), &vec![0]);

        assert!(poly1.clone().mul(&poly2).is_equal(&expected));
    }

    #[test]
    fn should_eval_multipolynomial() {
        let p = BigIntElement(BigInt::from(101));
        let g = BigIntElement(BigInt::from(0));
        let f = Rc::new(Field::new(p, g));

        // 4x + 2y^2 + 9
        let mut poly = MPolynomial::new(&f);
        poly.term(&BigIntElement(BigInt::from(4)), &vec![1]);
        poly.term(&BigIntElement(BigInt::from(2)), &vec![0, 2]);
        poly.term(&BigIntElement(BigInt::from(9)), &vec![0]);

        // x: 50, y: 20
        // 4*50 + 2*20^2 + 9
        let expected = f.modd(BigIntElement(BigInt::from(1009)));

        assert_eq!(
            poly.eval(&[
                BigIntElement(BigInt::from(50)),
                BigIntElement(BigInt::from(20))
            ]),
            expected
        );
    }

    #[test]
    fn should_eval_symbolic_multipolynomial() {
        let p = BigIntElement(BigInt::from(3221225473_u32));
        let g = BigIntElement(BigInt::from(5));
        let f = Rc::new(Field::new(p, g));

        // 4x + 2y^2 + 9
        let mut poly = MPolynomial::new(&f);
        poly.term(&BigIntElement(BigInt::from(4)), &vec![1]);
        poly.term(&BigIntElement(BigInt::from(2)), &vec![0, 2]);
        poly.term(&BigIntElement(BigInt::from(9)), &vec![0]);

        let mut x = Polynomial::new(&f);
        x.term(&BigIntElement(BigInt::from(5)), 2);

        let mut y = Polynomial::new(&f);
        y.term(&BigIntElement(BigInt::from(9)), 6);

        // x: 5x^2, y: 9x^6
        // 20x^2 + 162x^12 + 9
        let mut expected = Polynomial::new(&f);
        expected.term(&BigIntElement(BigInt::from(20)), 2);
        expected.term(&BigIntElement(BigInt::from(162)), 12);
        expected.term(&BigIntElement(BigInt::from(9)), 0);

        assert!(poly.eval_symbolic(&vec!(x, y)).is_equal(&expected));
    }

    #[test]
    fn should_make_multipolynomial_from_polynomial() {
        let p = BigIntElement(BigInt::from(101));
        let g = BigIntElement(BigInt::from(0));
        let f = Rc::new(Field::new(p, g));

        // 4x + 8x^2 + 9
        let mut poly = Polynomial::new(&f);
        poly.term(&BigIntElement(BigInt::from(4)), 1);
        poly.term(&BigIntElement(BigInt::from(8)), 2);
        poly.term(&BigIntElement(BigInt::from(9)), 0);

        let mpoly = MPolynomial::from_poly(&poly);

        for i in 0..20 {
            let v = BigIntElement(BigInt::from(i));
            assert_eq!(mpoly.eval(&vec!(v.clone())), poly.eval(&v));
        }
    }
}
