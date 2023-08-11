use num_bigint::BigInt;
use num_bigint::Sign;

#[derive(Clone)]
pub struct Field {
  p: BigInt,
  g: BigInt
}

impl Field {
  pub fn new(p: BigInt, g: BigInt) -> Field {
    Field {
      p, g
    }
  }

  pub fn bigintf(val: i32) -> BigInt {
    if val < 0 {
      BigInt::new(Sign::Minus, vec!(val.abs().try_into().unwrap()))
    } else {
      BigInt::new(Sign::Plus, vec!(val.try_into().unwrap()))
    }
  }

  pub fn biguintf(val: u32) -> BigInt {
    BigInt::new(Sign::Plus, vec!(val))
  }

  pub fn bigint(&self, val: i32) -> BigInt {
    self.modd(Field::bigintf(val))
  }

  pub fn biguint(&self, val: u32) -> BigInt {
    self.modd(BigInt::new(Sign::Plus, vec!(val)))
  }

  pub fn zero() -> BigInt {
    BigInt::new(Sign::Plus, vec!(0))
  }

  pub fn one() -> BigInt {
    BigInt::new(Sign::Plus, vec!(1))
  }

  pub fn modd(&self, v: BigInt) -> BigInt {
    if v < Field::zero() {
      let mut v_m = v.clone();
      while v_m < Field::zero() {
        v_m += &self.p;
      }
      return v_m;
    } else if v > self.p {
      return v % &self.p;
    }
    v
  }

  pub fn add(&self, v1: &BigInt, v2: &BigInt) -> BigInt {
    self.modd(v1 + v2)
  }

  pub fn mul(&self, v1: &BigInt, v2: &BigInt) -> BigInt {
    self.modd(v1 * v2)
  }

  pub fn sub(&self, v1: &BigInt, v2: &BigInt) -> BigInt {
    self.modd(v1 - v2)
  }

  pub fn neg(&self, v: &BigInt) -> BigInt {
    self.mul(v, &(&self.p-&Field::one()))
  }

  // exponent should always be > 0
  pub fn exp(&self, v: &BigInt, e: &BigInt) -> BigInt {
    self.modd(v.clone()).modpow(e, &self.p)
  }

  pub fn generator(&self, size: &BigInt) -> BigInt {
    if size >= &self.p {
      panic!("requested subgroup is larger than field");
    }
    let numer = &self.p - Field::one();
    let exp = &numer / size;
    if exp.clone() * size != numer {
      panic!("subgroup is not a divisor of field");
    }
    self.exp(&self.g, &exp)
  }

  pub fn inv(&self, v: &BigInt) -> BigInt {
    let mut val = self.modd(v.clone());
    if val == Field::zero() {
      panic!("divide by zero");
    }
    let one = Field::one();
    let mut y = Field::zero();
    let mut x = Field::one();
    let mut f = self.p.clone();

    while val > one {
      let q = val.clone() / f.clone();
      let mut t = f.clone();
      f = val % f;
      val = t;
      t = y.clone();
      y = x - q * y;
      x = t;
    }
    self.modd(x)
  }

}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn should_make_bigint() {
    let p = Field::bigintf(101);
    let f = Field::new(p, Field::bigintf(0));
    assert_eq!(f.bigint(0), BigInt::new(Sign::Minus, vec!(0)));
    assert_eq!(f.bigint(-29), BigInt::new(Sign::Plus , vec!(72)));
    assert_eq!(f.bigint(32), BigInt::new(Sign::Plus, vec!(32)));
    assert_eq!(f.biguint(0), BigInt::new(Sign::Minus, vec!(0)));
    assert_eq!(f.biguint(32), BigInt::new(Sign::Plus, vec!(32)));
  }

  #[test]
  fn should_add_two_elements() {
    let p = Field::bigintf(101);
    let g = Field::bigintf(0);
    let f = Field::new(p, g);

    let x = f.bigint(40);
    let y = f.bigint(90);

    assert_eq!(f.add(&x, &y), f.bigint(29));
  }

  #[test]
  fn should_mul_two_elements() {
    let p = Field::bigintf(101);
    let g = Field::bigintf(0);
    let f = Field::new(p, g);

    let x = f.bigint(40);
    let y = f.bigint(90);

    assert_eq!(f.mul(&x, &y), f.bigint(65));
  }

  #[test]
  fn should_sub_two_elements() {
    let p = Field::bigintf(101);
    let g = Field::bigintf(0);
    let f = Field::new(p, g);

    let x = f.bigint(2);
    let y = f.bigint(20);

    assert_eq!(f.sub(&x, &y), f.biguint(83));
  }

  #[test]
  fn should_get_generator() {
    let p = Field::biguintf(3221225473);
    let f_g = Field::bigintf(5);
    let f = Field::new(p, f_g);

    for i in 1..10 {
      let g = f.generator(&f.biguint(u32::pow(2, i)));
      assert_eq!(f.exp(&g, &f.biguint(u32::pow(2, i))), f.bigint(1));
    }
  }

  #[test]
  fn should_get_inverse() {
    let p = Field::biguintf(3221225473);
    let f_g = Field::bigintf(5);
    let f = Field::new(p, f_g);

    for i in 1..99 {
      let v = f.biguint(i);
      let inv = f.inv(&v);
      assert_eq!(f.mul(&inv, &v), Field::one());
    }
  }
}
