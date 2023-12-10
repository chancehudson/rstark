use crate::{field::Field, FieldElement};
use std::rc::Rc;

#[derive(Clone, Debug)]
pub struct Polynomial<T: FieldElement> {
    field: Rc<Field<T>>,
    coefs: Vec<T>,
}

impl<T: FieldElement> Polynomial<T> {
    pub fn field(&self) -> &Rc<Field<T>> {
        &self.field
    }

    pub fn new(f: &Rc<Field<T>>) -> Polynomial<T> {
        Polynomial {
            field: Rc::clone(f),
            coefs: Vec::new(),
        }
    }

    pub fn coefs(&self) -> &Vec<T> {
        &self.coefs
    }

    pub fn degree(&self) -> usize {
        let zero = self.field().zero();
        for i in 0..self.coefs.len() {
            let index = (self.coefs.len() - 1) - i;
            if self.coefs[index] != zero {
                return index;
            }
        }
        0
    }

    pub fn is_zero(&self) -> bool {
        let zero = self.field().zero();
        for i in 0..self.coefs.len() {
            if self.coefs[i] != zero {
                return false;
            }
        }
        true
    }

    pub fn is_equal(&self, poly: &Polynomial<T>) -> bool {
        let zero = self.field().zero();
        let larger;
        let smaller;
        if poly.coefs().len() > self.coefs().len() {
            larger = poly.coefs();
            smaller = self.coefs();
        } else {
            smaller = poly.coefs();
            larger = self.coefs();
        }
        for i in 0..smaller.len() {
            if larger.get(i) != smaller.get(i) {
                return false;
            }
        }
        for i in smaller.len()..larger.len() {
            if larger.get(i) != Some(&zero) {
                return false;
            }
        }
        true
    }

    pub fn term(&mut self, coef: &T, exp: u32) -> &Self {
        let s = usize::try_from(exp).unwrap();
        // expand to one longer to handle 0 exponents
        if self.coefs.len() < s + 1 {
            self.coefs.resize(s + 1, self.field().zero());
        }
        self.coefs[s] = self.field.add(&self.coefs[s], coef);
        self
    }

    pub fn add(&mut self, poly: &Polynomial<T>) -> &Self {
        for i in 0..self.coefs().len() {
            if i >= poly.coefs().len() {
                break;
            }
            self.coefs[i] = self.field.add(&self.coefs[i], &poly.coefs()[i]);
        }
        for i in self.coefs.len()..poly.coefs().len() {
            self.coefs.push(poly.coefs()[i].clone());
        }
        self
    }

    pub fn sub(&mut self, poly: &Polynomial<T>) -> &Self {
        for i in 0..self.coefs().len() {
            if i >= poly.coefs().len() {
                break;
            }
            self.coefs[i] = self.field.sub(&self.coefs[i], &poly.coefs()[i]);
        }
        for i in self.coefs.len()..poly.coefs().len() {
            self.coefs.push(self.field.neg(&poly.coefs()[i]));
        }
        self
    }

    // a fast method for multiplying by a single term polynomial
    // with a coefficient of 1
    // e.g. multiplying by x^5
    pub fn shift_and_clone(&self, degree: u32) -> Self {
        let degree_usize = usize::try_from(degree).unwrap();
        let mut shifted_coefs = vec![self.field().zero(); degree_usize];
        shifted_coefs.extend(self.coefs.clone());
        let mut out = Polynomial::new(&self.field);
        out.coefs = shifted_coefs;
        out
    }

    pub fn mul(&mut self, poly: &Polynomial<T>) -> &Self {
        let mut out = Vec::new();
        out.resize(self.coefs.len() + poly.coefs().len(), self.field().zero());
        for i in 0..poly.coefs().len() {
            // self.mul_term(&poly.coefs()[i], i);
            for j in 0..self.coefs.len() {
                // combine the exponents
                let e = j + i;
                out[e] = self
                    .field
                    .add(&out[e], &self.field.mul(&self.coefs[j], &poly.coefs()[i]));
            }
        }
        self.coefs = out;
        self.trim();
        self
    }

    pub fn mul_scalar(&mut self, v: &T) -> &Self {
        for i in 0..self.coefs.len() {
            self.coefs[i] = self.field.mul(v, &self.coefs[i]);
        }
        self
    }

    pub fn exp(&mut self, v: usize) -> &Self {
        let mut out = Polynomial::new(&self.field);
        out.term(&self.field().one(), 0);
        for _ in 0..v {
            out.mul(self);
        }
        self.coefs = out.coefs;
        self
    }

    // if we're scaling the polynomial using a generator point or similar
    // we probably already have a list of the exponents laying around
    pub fn scale_precalc(&mut self, _v: &T, exps: &Vec<T>) -> &Self {
        self.coefs = self
            .coefs
            .iter()
            .enumerate()
            .map(|(exp, coef)| self.field.mul(&exps[exp], coef))
            .collect();
        self
    }

    pub fn scale(&mut self, v: T) -> &Self {
        self.coefs = self
            .coefs
            .iter()
            .enumerate()
            .map(|(exp, coef)| {
                self.field.mul(
                    &self
                        .field
                        .exp(&v, &T::from_u32(exp as u32, self.field().p())),
                    coef,
                )
            })
            .collect();
        self
    }

    // compose `poly` into `this`
    pub fn compose(&mut self, poly: &Polynomial<T>) -> &Self {
        let mut out = Polynomial::new(&self.field);
        for (exp, coef) in self.coefs.iter().enumerate() {
            let mut p = poly.clone();
            p.exp(exp);
            p.mul_scalar(coef);
            out.add(&p);
        }
        self.coefs = out.coefs;
        self
    }

    // trim trailing zero coefficient
    pub fn trim(&mut self) {
        let mut new_len = self.coefs.len();
        let zero = self.field().zero();
        for i in (0..self.coefs.len()).rev() {
            if self.coefs[i] != zero {
                break;
            }
            new_len = i;
        }
        self.coefs.resize(new_len, zero);
    }

    // using horners method
    // https://en.wikipedia.org/wiki/Horner%27s_method
    pub fn eval(&self, v: &T) -> T {
        if self.coefs.is_empty() {
            return self.field().zero();
        }
        if self.coefs.len() == 1 {
            return self.coefs[0].clone();
        }
        let mut out = self.field.mul(v, &self.coefs[self.coefs.len() - 1]);
        for coef in self.coefs[1..(self.coefs.len() - 1)].iter().rev() {
            out = self.field.mul(v, &self.field.add(coef, &out));
        }
        out = self.field.add(&out, &self.coefs[0]);
        out
    }

    pub fn eval_batch(&self, vals: &Vec<T>) -> Vec<T> {
        vals.iter().map(|v| self.eval(v)).collect()
    }

    // https://en.wikipedia.org/wiki/Polynomial_evaluation#Multipoint_evaluation
    pub fn eval_batch_fast(&self, vals: &Vec<T>) -> Vec<T> {
        self.eval_batch_fast_slice(&vals[0..])
    }

    pub fn eval_batch_fast_slice(&self, vals: &[T]) -> Vec<T> {
        if vals.len() < 8 {
            return self.eval_batch(&vals.to_vec());
        }
        self.eval_batch_fast_(&vals[0..])
    }

    fn eval_batch_fast_(&self, vals: &[T]) -> Vec<T> {
        if vals.is_empty() {
            return Vec::new();
        }
        if vals.len() == 1 {
            return vec![self.eval(&vals[0])];
        }
        let half = vals.len() >> 1;

        let left_zeroifier = Polynomial::zeroifier_fft_slice(&vals[0..half], &self.field);
        let right_zeroifier = Polynomial::zeroifier_fft_slice(&vals[half..], &self.field);
        let (_, left_r) = self.div(&left_zeroifier);
        let (_, right_r) = self.div(&right_zeroifier);

        let mut left = left_r.eval_batch_fast_(&vals[0..half]);
        let right = right_r.eval_batch_fast_(&vals[half..]);
        left.extend(right);
        left
    }

    pub fn eval_batch_coset(&self, offset: &T, size: u32) -> Vec<T> {
        let mut scaled = self.clone();
        let offset_domain = self.field.domain(offset, size);
        scaled.scale_precalc(offset, &offset_domain);
        let (generator, _) = self.field.generator_cache(&size);
        let domain = self.field.domain(&generator, size);
        Self::eval_fft(scaled.coefs(), &domain, &self.field)
    }

    pub fn eval_batch_batch_coset(
        polys: Vec<Polynomial<T>>,
        offset: &T,
        size: u32,
        field: &Rc<Field<T>>,
    ) -> Vec<Vec<T>> {
        let offset_domain = field.domain(offset, size);
        let scaled = polys
            .iter()
            .map(|poly| {
                let mut p = poly.clone();
                p.scale_precalc(offset, &offset_domain);
                p
            })
            .collect::<Vec<Polynomial<T>>>();
        let (generator, _) = field.generator_cache(&size);
        let domain = field.domain(&generator, size);
        Self::eval_fft_batch(&scaled, &domain, field)
    }

    // Evaluate a polynomial over a multiplicative subgroup
    // domain cannot be a coset
    pub fn eval_batch_fft(&self, domain: &Vec<T>) -> Vec<T> {
        Self::eval_fft(self.coefs(), domain, &self.field)
    }

    pub fn eval_fft_batch(
        polys: &Vec<Polynomial<T>>,
        domain: &Vec<T>,
        field: &Rc<Field<T>>,
    ) -> Vec<Vec<T>> {
        let mut out = Vec::new();
        for _ in 0..polys.len() {
            let mut v = vec![field.zero(); domain.len()];
            v.reserve(5);
            out.push(v);
        }
        let mut scratch1 = field.zero();
        let mut scratch2 = field.zero();
        Self::eval_fft_batch_(
            polys,
            domain,
            field,
            1,
            0,
            0,
            domain.len() / 2,
            &mut out,
            &mut scratch1,
            &mut scratch2,
        );
        out
    }

    pub fn eval_fft(coefs: &Vec<T>, domain: &Vec<T>, field: &Rc<Field<T>>) -> Vec<T> {
        let mut out = vec![T::from_params(field.p()); domain.len()];
        Self::eval_fft_(coefs, domain, field, 1, 0, 0, domain.len() / 2, &mut out);
        out
    }

    pub fn eval_fft_batch_(
        polys: &[Polynomial<T>],
        domain: &[T],
        field: &Rc<Field<T>>,
        slice_len: usize,
        offset: usize,
        left_dest: usize,
        right_dest: usize,
        out: &mut Vec<Vec<T>>,
        scratch1: &mut T,
        scratch2: &mut T,
    ) {
        let out_size = domain.len() / slice_len;

        if out_size == 1 {
            for (i, p) in polys.iter().enumerate() {
                if let Some(v) = p.coefs().get(offset) {
                    out[i][left_dest].clone_from(v);
                }
            }
            // otherwise the coef is 0 which is the default value in `out`
            return;
        }

        Self::eval_fft_batch_(
            polys,
            domain,
            field,
            slice_len * 2,
            offset,
            left_dest,
            left_dest + out_size / 4,
            out,
            scratch1,
            scratch2,
        );
        Self::eval_fft_batch_(
            polys,
            domain,
            field,
            slice_len * 2,
            offset + slice_len,
            right_dest,
            right_dest + out_size / 4,
            out,
            scratch1,
            scratch2,
        );
        // let mut scratch = BigInt::from(0);

        for i in 0..(out_size / 2) {
            for j in 0..out.len() {
                let left_index = left_dest + i;
                let right_index = right_dest + i;
                // use this complicated scratch system to avoid
                // heap (de)allocations
                *scratch1 = field.mul(&out[j][right_index], &domain[i * slice_len]);
                *scratch2 = out[j][left_index].sub(scratch1);
                out[j][left_index] = out[j][left_index].add(scratch1);
                out[j][right_index].clone_from(scratch2);
            }
        }
    }

    // fft implemented by mutating a single `out` vec
    // this saves the cost of alloc/freeing vectors
    // during recursion
    //
    // measured improvement of ~10% compared to previous
    // approach
    pub fn eval_fft_(
        coefs: &Vec<T>,
        domain: &Vec<T>,
        field: &Rc<Field<T>>,
        slice_len: usize,
        offset: usize,
        left_dest: usize,
        right_dest: usize,
        out: &mut Vec<T>,
    ) {
        let out_size = domain.len() / slice_len;

        if out_size == 1 {
            if let Some(v) = coefs.get(offset) {
                out[left_dest] = v.clone();
            }
            // otherwise the coef is 0 which is the default value in `out`
            return;
        }

        Self::eval_fft_(
            coefs,
            domain,
            field,
            slice_len * 2,
            offset,
            left_dest,
            left_dest + out_size / 4,
            out,
        );
        Self::eval_fft_(
            coefs,
            domain,
            field,
            slice_len * 2,
            offset + slice_len,
            right_dest,
            right_dest + out_size / 4,
            out,
        );

        for i in 0..(out_size / 2) {
            let left_out;
            let right_out;
            {
                let x = &out[left_dest + i];
                let y = &out[right_dest + i];
                let y_root = field.mul(y, &domain[i * slice_len]);
                left_out = x.add(&y_root);
                right_out = x.sub(&y_root);
            }
            out[left_dest + i] = left_out;
            out[right_dest + i] = right_out;
        }
    }

    pub fn eval_fft_inv(coefs: &Vec<T>, domain_inv: &Vec<T>, field: &Rc<Field<T>>) -> Vec<T> {
        if coefs.len() == 1 {
            return vec![coefs[0].clone()];
        }
        let len_inv = field.inv(&T::from_u32(u32::try_from(coefs.len()).unwrap(), field.p()));
        let out = Self::eval_fft(coefs, domain_inv, field);
        out.into_iter().map(|v| field.mul(&v, &len_inv)).collect()
    }

    pub fn mul_fft(
        poly1: &Polynomial<T>,
        poly2: &Polynomial<T>,
        field: &Rc<Field<T>>,
    ) -> Polynomial<T> {
        if poly1.degree() < 4 || poly2.degree() < 4 {
            let mut o = poly1.clone();
            o.mul(poly2);
            return o;
        }
        let out_degree = 2 * std::cmp::max(poly1.degree(), poly2.degree());
        let domain_size = 2_u32.pow(u32::try_from(out_degree).unwrap().ilog2() + 1);
        let (generator, generator_inv) = field.generator_cache(&domain_size);
        let domain = field.domain(&generator, domain_size);

        let x1 = Self::eval_fft(poly1.coefs(), &domain, field);
        let x2 = Self::eval_fft(poly2.coefs(), &domain, field);

        let mut x3 = Vec::new();
        for i in 0..domain.len() {
            x3.push(field.mul(&x1[i], &x2[i]));
        }
        let domain_inv = field.domain(&generator_inv, domain_size);

        let out = Self::eval_fft_inv(&x3, &domain_inv, field);
        let mut p = Polynomial {
            field: Rc::clone(field),
            coefs: out,
        };
        p.trim();
        p
    }

    pub fn div_coset(
        poly1: &Polynomial<T>,
        poly2: &Polynomial<T>,
        offset: &T,
        generator: &T,
        size: u32,
        field: &Rc<Field<T>>,
    ) -> Polynomial<T> {
        if poly1.is_zero() {
            return Polynomial::new(field);
        }
        if poly2.is_zero() {
            panic!("divide by 0");
        }
        if poly2.degree() > poly1.degree() {
            panic!("divisor cannot have degree larger than dividend");
        }
        let degree = poly1.degree();

        let mut g = generator.clone();
        let mut order = size;

        while u32::try_from(degree).unwrap() < order >> 1 {
            g = field.mul(&g, &g);
            order >>= 1;
        }

        let g_inv = field.inv(&g);
        let domain = field.domain(&g, order);
        let domain_inv = field.domain(&g_inv, order);
        let offset_domain = field.domain(offset, order);

        let mut poly1_scaled = poly1.clone();
        poly1_scaled.scale_precalc(offset, &offset_domain);
        let mut poly2_scaled = poly2.clone();
        poly2_scaled.scale_precalc(offset, &offset_domain);

        let poly1_codeword = Self::eval_fft(poly1_scaled.coefs(), &domain, field);
        let poly2_codeword = Self::eval_fft(poly2_scaled.coefs(), &domain, field);

        let poly2_codeword_inv = field.inv_batch(&poly2_codeword);

        let out = poly1_codeword
            .iter()
            .enumerate()
            .map(|(i, val)| field.mul(val, &poly2_codeword_inv[i]))
            .collect();

        let scaled_coefs = Self::eval_fft_inv(&out, &domain_inv, field);
        let mut scaled_poly = Polynomial::new(field);
        for i in 0..(poly1.degree() - poly2.degree() + 1) {
            scaled_poly.term(&scaled_coefs[i], u32::try_from(i).unwrap());
        }
        let offset_inv = field.inv(offset);
        let offset_domain_inv = field.domain(&offset_inv, order);
        scaled_poly.scale_precalc(&offset_inv, &offset_domain_inv);
        scaled_poly
    }

    // remove and return the largest non-zero coefficient
    // coef, exp
    pub fn pop_term(&mut self) -> (T, usize) {
        let zero = self.field().zero();
        for i in 0..self.coefs.len() {
            let index = (self.coefs.len() - 1) - i;
            if self.coefs[index] != zero {
                let out = self.coefs[index].clone();
                self.coefs[index] = zero;
                return (out, index);
            }
        }
        (zero, 0)
    }

    pub fn safe_div(&self, divisor: &Polynomial<T>) -> Polynomial<T> {
        let (q, r) = self.div(divisor);
        if !r.is_zero() {
            panic!("non-zero remainder in division");
        }
        q
    }

    pub fn div(&self, divisor: &Polynomial<T>) -> (Polynomial<T>, Polynomial<T>) {
        if divisor.is_zero() {
            panic!("divide by zero");
        }
        let mut dclone = divisor.clone();
        let mut q = Polynomial::new(&self.field);
        let divisor_term = dclone.pop_term();
        let divisor_term_inv = self.field.inv(&divisor_term.0);
        let mut inter = self.clone();
        while inter.degree() >= divisor.degree() {
            let largest_term = inter.clone().pop_term();
            let new_coef = self.field.mul(&largest_term.0, &divisor_term_inv);
            let new_exp = largest_term.1 - divisor_term.1;
            q.term(&new_coef, new_exp.try_into().unwrap());
            let mut t = divisor.shift_and_clone(u32::try_from(new_exp).unwrap());
            t.mul_scalar(&new_coef);
            inter.sub(&t);
        }
        (q, inter)
    }

    pub fn lagrange(x_vals: &Vec<T>, y_vals: &Vec<T>, field: &Rc<Field<T>>) -> Polynomial<T> {
        if x_vals.len() != y_vals.len() {
            panic!("lagrange mismatch x/y array length");
        }
        let mut numerator = Polynomial::new(field);
        numerator.term(&field.one(), 0);
        for v in x_vals {
            let mut poly = Polynomial::new(field);
            poly.term(&field.bigint(1), 1);
            poly.term(&field.neg(v), 0);
            numerator.mul(&poly);
        }
        let mut polynomials = Vec::new();
        for i in 0..x_vals.len() {
            let mut denominator = field.one();
            for j in 0..x_vals.len() {
                if i == j {
                    continue;
                }
                denominator = field.mul(&denominator, &x_vals[i].sub(&x_vals[j]));
            }
            let mut n = Polynomial::new(field);
            n.term(&field.bigint(1), 1);
            n.term(&field.neg(&x_vals[i]), 0);
            n.mul_scalar(&denominator);
            let poly = numerator.safe_div(&n);
            polynomials.push(poly);
        }
        let mut out = Polynomial::new(field);
        for i in 0..x_vals.len() {
            out.add(polynomials[i].mul_scalar(&y_vals[i]));
        }
        out
    }

    pub fn interpolate_fft(
        x_vals: &Vec<T>,
        y_vals: &Vec<T>,
        field: &Rc<Field<T>>,
    ) -> Polynomial<T> {
        Self::interpolate_fft_slice(&x_vals[0..], &y_vals[0..], field)
    }

    pub fn interpolate_fft_slice(
        x_vals: &[T],
        y_vals: &[T],
        field: &Rc<Field<T>>,
    ) -> Polynomial<T> {
        if x_vals.len() != y_vals.len() {
            panic!("x/y len mismatch");
        }
        if x_vals.is_empty() {
            return Polynomial::new(field);
        }
        if x_vals.len() == 1 {
            let mut p = Polynomial::new(field);
            p.term(&y_vals[0], 0);
            return p;
        }
        let half = x_vals.len() >> 1;

        let left_zeroifier = Self::zeroifier_fft_slice(&x_vals[0..half], field);
        let right_zeroifier = Self::zeroifier_fft_slice(&x_vals[half..], field);

        let left_offset = right_zeroifier.eval_batch_fast_slice(&x_vals[0..half]);
        let right_offset = left_zeroifier.eval_batch_fast_slice(&x_vals[half..]);

        let left_offset_inv = field.inv_batch(&left_offset);
        let right_offset_inv = field.inv_batch(&right_offset);

        let left_targets = y_vals[0..half]
            .iter()
            .enumerate()
            .map(|(i, v)| field.mul(v, &left_offset_inv[i]))
            .collect::<Vec<T>>();
        let right_targets = y_vals[half..]
            .iter()
            .enumerate()
            .map(|(i, v)| field.mul(v, &right_offset_inv[i]))
            .collect::<Vec<T>>();

        let mut left_interpolant =
            Self::interpolate_fft_slice(&x_vals[0..half], &left_targets[0..], field);
        let mut right_interpolant =
            Self::interpolate_fft_slice(&x_vals[half..], &right_targets[0..], field);

        left_interpolant.mul(&right_zeroifier);
        right_interpolant.mul(&left_zeroifier);
        left_interpolant.add(&right_interpolant);
        left_interpolant
    }

    pub fn interpolate_fft_batch(
        x_vals: &[T],
        y_vals: &[Vec<T>],
        field: &Rc<Field<T>>,
    ) -> Vec<Polynomial<T>> {
        if x_vals.is_empty() {
            return vec![Polynomial::new(field)];
        }
        if x_vals.len() == 1 {
            return y_vals
                .iter()
                .map(|vals| {
                    let mut p = Polynomial::new(field);
                    p.term(&vals[0], 0);
                    p
                })
                .collect();
        }
        let half = x_vals.len() >> 1;

        let left_zeroifier = Self::zeroifier_fft_slice(&x_vals[0..half], field);
        let right_zeroifier = Self::zeroifier_fft_slice(&x_vals[half..], field);

        let left_offset = right_zeroifier.eval_batch_fast_slice(&x_vals[0..half]);
        let right_offset = left_zeroifier.eval_batch_fast_slice(&x_vals[half..]);

        let left_offset_invs: Vec<T> = left_offset.iter().map(|v| field.inv(v)).collect();
        let right_offset_invs: Vec<T> = right_offset.iter().map(|v| field.inv(v)).collect();

        let left_targets: Vec<Vec<T>> = y_vals
            .iter()
            .map(|vals| {
                vals[0..half]
                    .iter()
                    .enumerate()
                    .map(|(i, v)| field.mul(v, &left_offset_invs[i]))
                    .collect::<Vec<T>>()
            })
            .collect();
        let right_targets: Vec<Vec<T>> = y_vals
            .iter()
            .map(|vals| {
                vals[half..]
                    .iter()
                    .enumerate()
                    .map(|(i, v)| field.mul(v, &right_offset_invs[i]))
                    .collect::<Vec<T>>()
            })
            .collect();

        let left_interpolant =
            Self::interpolate_fft_batch(&x_vals[0..half], &left_targets[0..], field);
        let right_interpolant =
            Self::interpolate_fft_batch(&x_vals[half..], &right_targets[0..], field);

        let mut out = Vec::new();
        for i in 0..left_interpolant.len() {
            let mut left = Polynomial::mul_fft(&left_interpolant[i], &right_zeroifier, field);
            let right = Polynomial::mul_fft(&right_interpolant[i], &left_zeroifier, field);
            left.add(&right);
            out.push(left);
        }
        out
    }

    pub fn test_colinearity_batch(
        x_vals_arr: &[Vec<T>],
        y_vals_arr: &[Vec<T>],
        field: &Rc<Field<T>>,
    ) -> bool {
        let mut to_inv = Vec::new();
        for x_vals in x_vals_arr {
            to_inv.push(field.sub(&x_vals[1], &x_vals[0]));
            to_inv.push(field.sub(&x_vals[2], &x_vals[1]));
        }
        let inverted = field.inv_batch(&to_inv);
        for (i, y_vals) in y_vals_arr.iter().enumerate() {
            let x_diff_inv_1 = &inverted[2 * i];
            let x_diff_inv_2 = &inverted[2 * i + 1];
            let y_diff_1 = field.sub(&y_vals[1], &y_vals[0]);
            let y_diff_2 = field.sub(&y_vals[2], &y_vals[1]);
            let slope_1 = field.mul(&y_diff_1, x_diff_inv_1);
            let slope_2 = field.mul(&y_diff_2, x_diff_inv_2);
            if slope_1 != slope_2 {
                return false;
            }
        }
        true
    }

    pub fn test_colinearity(x_vals: &Vec<T>, y_vals: &Vec<T>, field: &Rc<Field<T>>) -> bool {
        let x_diff_1 = field.sub(&x_vals[1], &x_vals[0]);
        let y_diff_1 = field.sub(&y_vals[1], &y_vals[0]);
        let slope_1 = field.div(&y_diff_1, &x_diff_1);
        let x_diff_2 = field.sub(&x_vals[2], &x_vals[1]);
        let y_diff_2 = field.sub(&y_vals[2], &y_vals[1]);
        let slope_2 = field.div(&y_diff_2, &x_diff_2);
        slope_1 == slope_2
    }

    pub fn zeroifier_fft(points: &Vec<T>, field: &Rc<Field<T>>) -> Polynomial<T> {
        Self::zeroifier_fft_slice(&points[0..], field)
    }

    pub fn zeroifier_fft_slice(points: &[T], field: &Rc<Field<T>>) -> Polynomial<T> {
        if points.is_empty() {
            return Polynomial::new(field);
        }
        if points.len() == 1 {
            let mut p = Polynomial::new(field);
            p.term(&field.neg(&points[0]), 0);
            p.term(&field.one(), 1);
            return p;
        }
        let half = points.len() >> 1;
        let left = Self::zeroifier_fft_slice(&points[0..(half)], field);
        let right = Self::zeroifier_fft_slice(&points[(half)..], field);
        Self::mul_fft(&left, &right, field)
    }

    pub fn zeroifier(points: &Vec<T>, field: &Rc<Field<T>>) -> Polynomial<T> {
        Self::zeroifier_slice(&points[0..], field)
    }

    pub fn zeroifier_slice(points: &[T], field: &Rc<Field<T>>) -> Polynomial<T> {
        let mut out = Polynomial::new(field);
        out.term(&field.one(), 0);
        let mut x = Polynomial::new(field);
        x.term(&field.one(), 1);
        for p in points {
            out.mul(x.clone().term(&field.neg(p), 0));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;

    use crate::{to_crypto_element, to_crypto_params, BigIntElement};

    use super::*;

    #[test]
    fn should_test_colinearity() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        let mut poly = Polynomial::new(&f);
        poly.term(&f.bigint(-9), 0);
        poly.term(&f.bigint(2), 1);
        let mut x_vals = Vec::new();
        let mut y_vals = Vec::new();
        for i in 0..3 {
            x_vals.push(f.bigint(i));
            y_vals.push(poly.eval(&f.bigint(i)));
        }
        assert!(Polynomial::test_colinearity(&x_vals, &y_vals, &f));
        assert!(Polynomial::test_colinearity_batch(
            &vec!(x_vals),
            &vec!(y_vals),
            &f
        ));
    }

    #[test]
    fn should_compose_polynomial() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        let mut root = Polynomial::new(&f);
        root.term(&f.bigint(99), 0);
        root.term(&f.bigint(2), 1);
        root.term(&f.bigint(4), 2);

        let mut inpoly = Polynomial::new(&f);
        inpoly.term(&f.bigint(2), 2);
        inpoly.term(&f.bigint(12), 0);

        let mut expected = Polynomial::new(&f);
        expected.term(&f.bigint(99), 0);
        expected.add(&inpoly.clone().mul_scalar(&f.bigint(2)));
        {
            let mut i = inpoly.clone();
            i.exp(2);
            i.mul_scalar(&f.bigint(4));
            expected.add(&i);
        }

        assert!(root.compose(&inpoly).is_equal(&expected));
    }

    #[test]
    fn should_exp_polynomial() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        let mut poly = Polynomial::new(&f);
        poly.term(&f.bigint(2), 0);

        for i in 0..10 {
            let mut expected = Polynomial::new(&f);
            expected.term(&f.exp(&f.bigint(2), &f.bigint(i)), 0);
            assert!(poly.clone().exp(i as usize).is_equal(&expected));
        }
    }

    #[test]
    fn should_scale_polynomial() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        let mut poly = Polynomial::new(&f);
        poly.term(&f.bigint(1), 0);
        poly.term(&f.bigint(1), 1);
        poly.term(&f.bigint(1), 2);
        poly.term(&f.bigint(1), 3);
        poly.term(&f.bigint(1), 4);

        let mut expected = Polynomial::new(&f);
        expected.term(&f.bigint(1), 0);
        expected.term(&f.bigint(2), 1);
        expected.term(&f.bigint(4), 2);
        expected.term(&f.bigint(8), 3);
        expected.term(&f.bigint(16), 4);

        assert!(expected.is_equal(poly.scale(f.bigint(2))));
    }

    #[test]
    fn should_interpolate_lagrange() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        let size = 32;
        let gen = f.generator(f.bigint(size));
        let mut x_vals = Vec::new();
        let mut y_vals = Vec::new();
        for i in 0..size {
            x_vals.push(f.exp(&gen, &f.bigint(i)));
            y_vals.push(f.random());
        }

        let p = Polynomial::lagrange(&x_vals, &y_vals, &f);
        for i in 0..size {
            let index = i as usize;
            assert_eq!(p.eval(&x_vals[index]), y_vals[index]);
        }
    }

    #[test]
    #[should_panic]
    fn should_fail_to_div_by_zero() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        let mut poly = Polynomial::new(&f);
        poly.term(&f.bigint(1), 0);

        let mut divisor = Polynomial::new(&f);
        divisor.term(&f.bigint(0), 3);
        poly.div(&divisor);
    }

    #[test]
    fn should_divide_polynomial() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        let mut poly1 = Polynomial::new(&f);
        poly1.term(&f.bigint(1), 2);
        poly1.term(&f.bigint(2), 1);
        poly1.term(&f.bigint(-7), 0);

        let mut poly2 = Polynomial::new(&f);
        poly2.term(&f.bigint(1), 1);
        poly2.term(&f.bigint(-2), 0);

        let mut expected_q = Polynomial::new(&f);
        expected_q.term(&f.bigint(1), 1);
        expected_q.term(&f.bigint(4), 0);

        let mut expected_r = Polynomial::new(&f);
        expected_r.term(&f.bigint(1), 0);

        let (q, r) = poly1.div(&poly2);

        assert!(q.is_equal(&expected_q));
        assert!(r.is_equal(&expected_r));
    }

    #[test]
    fn should_eval_polynomial() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        // 9x^3 - 4x^2 - 20
        let mut poly = Polynomial::new(&f);
        poly.term(&f.bigint(9), 3);
        poly.term(&f.bigint(-4), 2);
        poly.term(&f.bigint(-20), 0);

        assert_eq!(f.bigint(-20), poly.eval(&f.bigint(0)));
        assert_eq!(f.bigint(-15), poly.eval(&f.bigint(1)));
    }

    #[test]
    fn should_eval_polynomial_with_batch_fast() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        let mut poly = Polynomial::new(&f);
        for i in 0..50 {
            poly.term(&f.random(), i);
        }

        let size = 2_u32.pow(7);
        let mut g_domain = Vec::new();
        for i in 0..size {
            g_domain.push(f.biguint(i));
        }

        let actual = poly.eval_batch(&g_domain);
        let out = poly.eval_batch_fast(&g_domain);
        for i in 0..usize::try_from(size).unwrap() {
            assert_eq!(actual[i], out[i]);
        }
    }

    #[test]
    fn should_eval_polynomial_with_fft() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        let mut poly = Polynomial::new(&f);
        for i in 0..50 {
            poly.term(&f.random(), i);
        }

        let size = 2_u32.pow(8);
        let sub_g = f.generator(f.biguint(size));
        let mut g_domain = Vec::new();
        for i in 0..size {
            g_domain.push(f.exp(&sub_g, &f.biguint(i)));
        }

        let actual = poly.eval_batch(&g_domain);
        let out = poly.eval_batch_fft(&g_domain);
        for i in 0..usize::try_from(size).unwrap() {
            assert_eq!(actual[i], out[i]);
        }
    }

    #[test]
    fn should_check_polynomial_equality() {
        let p = to_crypto_params(BigIntElement(BigInt::from(101)));
        let g = to_crypto_element(BigIntElement(BigInt::from(0)), &p);
        let f = Rc::new(Field::new(g));

        let mut poly1 = Polynomial::new(&f);
        poly1.term(&f.bigint(-20), 0);
        poly1.term(&f.bigint(2), 2);
        let mut poly2 = Polynomial::new(&f);
        poly2.term(&f.bigint(81), 0);
        poly2.term(&f.bigint(2), 2);
        poly2.term(&f.bigint(0), 5);
        assert!(poly1.is_equal(&poly2));
    }

    #[test]
    fn should_add_polynomials() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        // 2x^2 - 20
        let mut poly1 = Polynomial::new(&f);
        poly1.term(&f.bigint(-20), 0);
        poly1.term(&f.bigint(2), 2);

        // 9x^3 - 4x^2 - 20
        let mut poly2 = Polynomial::new(&f);
        poly2.term(&f.bigint(9), 3);
        poly2.term(&f.bigint(-4), 2);
        poly2.term(&f.bigint(-20), 0);

        // expect 9x^3 - 2x^2 - 40
        let mut out1 = poly1.clone();
        out1.add(&poly2);
        let mut out2 = poly2.clone();
        out2.add(&poly1);

        let mut expected = Polynomial::new(&f);
        expected.term(&f.bigint(9), 3);
        expected.term(&f.bigint(-2), 2);
        expected.term(&f.bigint(-40), 0);
        assert!(expected.is_equal(&out1));
        assert!(expected.is_equal(&out2));
    }

    #[test]
    fn should_subtract_polynomials() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        // 2x^2 - 20
        let mut poly1 = Polynomial::new(&f);
        poly1.term(&f.bigint(-20), 0);
        poly1.term(&f.bigint(2), 2);

        // 9x^3 - 4x^2 - 20
        let mut poly2 = Polynomial::new(&f);
        poly2.term(&f.bigint(9), 3);
        poly2.term(&f.bigint(-4), 2);
        poly2.term(&f.bigint(-20), 0);

        let mut out1 = poly1.clone();
        out1.sub(&poly2);
        let mut out2 = poly2.clone();
        out2.sub(&poly1);

        let mut expected1 = Polynomial::new(&f);
        expected1.term(&f.bigint(-9), 3);
        expected1.term(&f.bigint(6), 2);
        assert!(expected1.is_equal(&out1));

        let mut expected2 = Polynomial::new(&f);
        expected2.term(&f.bigint(9), 3);
        expected2.term(&f.bigint(-6), 2);
        assert!(expected2.is_equal(&out2));
    }

    #[test]
    fn should_multiply_polynomials() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        // 2x^2 - 20
        let mut poly1 = Polynomial::new(&f);
        poly1.term(&f.bigint(-20), 0);
        poly1.term(&f.bigint(2), 2);

        // 9x^3 - 4x^2 - 20
        let mut poly2 = Polynomial::new(&f);
        poly2.term(&f.bigint(9), 3);
        poly2.term(&f.bigint(-4), 2);
        poly2.term(&f.bigint(-20), 0);

        let mut poly3 = poly1.clone();
        poly3.mul(&poly2);

        let mut expected = Polynomial::new(&f);
        expected.term(&f.bigint(18), 5);
        expected.term(&f.bigint(-8), 4);
        expected.term(&f.bigint(-180), 3);
        expected.term(&f.bigint(40), 2);
        expected.term(&f.bigint(400), 0);
        assert!(poly3.is_equal(&expected));
    }

    #[test]
    fn should_multiply_polynomials_fft() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        let mut poly1 = Polynomial::new(&f);
        let mut poly2 = Polynomial::new(&f);
        for i in 0..2 {
            poly1.term(&f.random(), i);
            poly2.term(&f.random(), i);
        }

        let mut expected = poly1.clone();
        expected.mul(&poly2);

        let out = Polynomial::mul_fft(&poly1, &poly2, &f);
        assert!(out.is_equal(&expected));
    }

    #[test]
    fn should_build_zeroifier_polynomial() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));

        let s = 128;
        let size = f.biguint(s);
        let domain_g = f.generator(size);
        let mut domain = Vec::new();
        for i in 0..s {
            domain.push(f.exp(&domain_g, &f.biguint(i)));
        }

        let poly = Polynomial::zeroifier(&domain, &f);
        for i in 0..s {
            assert_eq!(poly.eval(&domain[i as usize]), f.zero());
        }

        let mut is_zero = true;
        for i in poly.coefs() {
            if i != &f.zero() {
                is_zero = false;
            }
        }
        assert!(!is_zero);
    }

    #[test]
    fn should_build_zeroifier_polynomial_domain() {
        let p = to_crypto_params(BigIntElement(BigInt::from(3221225473_u32)));
        let g = to_crypto_element(BigIntElement(BigInt::from(5)), &p);
        let f = Rc::new(Field::new(g));
        let size = 128_u32;
        let generator = f.generator(f.biguint(size));
        let domain = f.domain(&generator, size);
        let zeroifier = Polynomial::zeroifier(&domain, &f);
        let mut zeroifier_fft = Polynomial::zeroifier_fft(&domain, &f);
        zeroifier_fft.trim();

        for v in domain {
            assert_eq!(f.zero(), zeroifier_fft.eval(&v));
            assert_eq!(f.zero(), zeroifier.eval(&v));
        }

        assert!(zeroifier.is_equal(&zeroifier_fft));
    }
}
