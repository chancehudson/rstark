use num_bigint::BigInt;
use crate::field::Field;
use std::rc::Rc;

#[derive(Clone)]
pub struct Polynomial {
  field: Rc<Field>,
  coefs: Vec<BigInt>
}

impl Polynomial {
  pub fn new(f: &Rc<Field>) -> Polynomial {
    Polynomial {
      field: Rc::clone(f),
      coefs: Vec::new()
    }
  }

  pub fn coefs(&self) -> &Vec<BigInt> {
    &self.coefs
  }

  pub fn degree(&self) -> usize {
    let zero = Field::zero();
    for i in 0..self.coefs.len() {
      let index = (self.coefs.len() - 1) - i;
      if self.coefs[index] != zero {
        return index;
      }
    }
    return 0;
  }

  pub fn is_zero(&self) -> bool {
    let zero = Field::zero();
    for i in 0..self.coefs.len() {
      if self.coefs[i] != zero {
        return false;
      }
    }
    true
  }

  pub fn is_equal(&self, poly: &Polynomial) -> bool {
    let zero = Field::zero();
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

  pub fn term(& mut self, coef: &BigInt, exp: u32) -> &Self {
    let s = usize::try_from(exp).unwrap();
    // expand to one longer to handle 0 exponents
    if self.coefs.len() < s+1 {
      self.coefs.resize(s+1, Field::zero());
    }
    self.coefs[s] = self.field.add(&self.coefs[s], &coef);
    self
  }

  pub fn add(& mut self, poly: &Polynomial) -> &Self {
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

  pub fn sub(& mut self, poly: &Polynomial) -> &Self {
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

  pub fn mul(& mut self, poly: &Polynomial) -> &Self {
    let mut out: Vec<BigInt> = Vec::new();
    out.resize(self.coefs.len() + poly.coefs().len(), Field::zero());
    for i in 0..poly.coefs().len() {
      // self.mul_term(&poly.coefs()[i], i);
      for j in 0..self.coefs.len() {
        // combine the exponents
        let e = j + i;
        out[e] = self.field.add(&out[e], &self.field.mul(&self.coefs[j], &poly.coefs()[i]));
      }
    }
    self.coefs = out;
    self.trim();
    self
  }

  // trim trailing zero coefficient
  pub fn trim(& mut self) {
    let mut new_len = self.coefs.len();
    let zero = Field::zero();
    for i in self.coefs.len()..0 {
      if self.coefs[i] != zero {
        break;
      }
      new_len = i;
    }
    self.coefs.resize(new_len, zero);
  }

  pub fn eval(&self, v: &BigInt) -> BigInt {
    let mut out = Field::zero();
    for i in 0..self.coefs.len() {
      out += &self.coefs[i] * self.field.exp(v, &self.field.bigint(i.try_into().unwrap()));
    }
    self.field.modd(out)
  }

  // remove and return the largest non-zero coefficient
  // coef, exp
  pub fn pop_term(& mut self) -> (BigInt, usize) {
    let zero = Field::zero();
    for i in self.coefs.len()..0 {
    println!("{}", i);
    }
    for i in 0..self.coefs.len() {
      let index = (self.coefs.len() - 1) - i;
      if self.coefs[index] != zero {
        let out = self.coefs[index].clone();
        self.coefs[index] = zero;
        return (out, index)
      }
    }
    return (zero, 0);
  }

  pub fn div(&self, divisor: &Polynomial) -> (Polynomial, Polynomial) {
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
      let mut t = Polynomial::new(&self.field);
      t.term(&new_coef, new_exp.try_into().unwrap());
      t.mul(divisor);
      inter.sub(&t);
    }
    (q, inter)
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  #[should_panic]
  fn should_fail_to_div_by_zero() {
    let p = Field::biguintf(3221225473);
    let g = Field::bigintf(5);
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    poly.term(&f.bigint(1), 0);

    let mut divisor = Polynomial::new(&f);
    divisor.term(&f.bigint(0), 3);
    poly.div(&divisor);
  }

  #[test]
  fn should_divide_polynomial() {
    let p = Field::biguintf(3221225473);
    let g = Field::bigintf(5);
    let f = Rc::new(Field::new(p, g));

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
    let p = Field::biguintf(3221225473);
    let g = Field::bigintf(5);
    let f = Rc::new(Field::new(p, g));

    // 9x^3 - 4x^2 - 20
    let mut poly = Polynomial::new(&f);
    poly.term(&f.bigint(9), 3);
    poly.term(&f.bigint(-4), 2);
    poly.term(&f.bigint(-20), 0);

    assert_eq!(f.bigint(-20), poly.eval(&f.bigint(0)));
    assert_eq!(f.bigint(-15), poly.eval(&f.bigint(1)));
  }

  #[test]
  fn should_check_polynomial_equality() {
    let p = Field::bigintf(101);
    let g = Field::bigintf(0);
    let f = Rc::new(Field::new(p, g));

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
    let p = Field::bigintf(101);
    let g = Field::bigintf(0);
    let f = Rc::new(Field::new(p, g));

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
    let p = Field::bigintf(101);
    let g = Field::bigintf(0);
    let f = Rc::new(Field::new(p, g));

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
    let p = Field::bigintf(101);
    let g = Field::bigintf(0);
    let f = Rc::new(Field::new(p, g));

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
}
