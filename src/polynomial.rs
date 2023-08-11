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

  pub fn term(& mut self, coef: BigInt, exp: u32) -> &Self {
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
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn should_check_polynomial_equality() {
    let p = Field::bigintf(101);
    let g = Field::bigintf(0);
    let f = Rc::new(Field::new(p, g));

    let mut poly1 = Polynomial::new(&f);
    poly1.term(f.bigint(-20), 0);
    poly1.term(f.bigint(2), 2);
    let mut poly2 = Polynomial::new(&f);
    poly2.term(f.bigint(81), 0);
    poly2.term(f.bigint(2), 2);
    poly2.term(f.bigint(0), 5);
    assert!(poly1.is_equal(&poly2));
  }

  #[test]
  fn should_add_polynomials() {
    let p = Field::bigintf(101);
    let g = Field::bigintf(0);
    let f = Rc::new(Field::new(p, g));

    // 2x^2 - 20
    let mut poly1 = Polynomial::new(&f);
    poly1.term(f.bigint(-20), 0);
    poly1.term(f.bigint(2), 2);

    // 9x^3 - 4x^2 - 20
    let mut poly2 = Polynomial::new(&f);
    poly2.term(f.bigint(9), 3);
    poly2.term(f.bigint(-4), 2);
    poly2.term(f.bigint(-20), 0);

    // expect 9x^3 - 2x^2 - 40
    let mut out1 = poly1.clone();
    out1.add(&poly2);
    let mut out2 = poly2.clone();
    out2.add(&poly1);

    let mut expected = Polynomial::new(&f);
    expected.term(f.bigint(9), 3);
    expected.term(f.bigint(-2), 2);
    expected.term(f.bigint(-40), 0);
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
    poly1.term(f.bigint(-20), 0);
    poly1.term(f.bigint(2), 2);

    // 9x^3 - 4x^2 - 20
    let mut poly2 = Polynomial::new(&f);
    poly2.term(f.bigint(9), 3);
    poly2.term(f.bigint(-4), 2);
    poly2.term(f.bigint(-20), 0);

    let mut out1 = poly1.clone();
    out1.sub(&poly2);
    let mut out2 = poly2.clone();
    out2.sub(&poly1);

    let mut expected1 = Polynomial::new(&f);
    expected1.term(f.bigint(-9), 3);
    expected1.term(f.bigint(6), 2);
    assert!(expected1.is_equal(&out1));

    let mut expected2 = Polynomial::new(&f);
    expected2.term(f.bigint(9), 3);
    expected2.term(f.bigint(-6), 2);
    assert!(expected2.is_equal(&out2));
  }
}
