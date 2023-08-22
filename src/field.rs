use num_bigint::{BigInt, Sign};
use rand::Rng;
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::sync::RwLock;

#[derive(Serialize, Deserialize)]
pub struct Field {
  p: BigInt,
  g: BigInt,
  // generator, size size keyed to contents
  group_cache: RwLock<HashMap<(BigInt, u32), Vec<BigInt>>>,
  // group size, offset keyed to contents
  coset_cache: RwLock<HashMap<(u32, BigInt), Vec<BigInt>>>
}

impl Field {
  pub fn new(p: BigInt, g: BigInt) -> Field {
    Field {
      p,
      g,
      group_cache: RwLock::new(HashMap::new()),
      coset_cache: RwLock::new(HashMap::new())
    }
  }

  pub fn bigint_to_u32(v: &BigInt) -> u32 {
    let (_, digits) = v.to_u32_digits();
    if digits.len() != 1 {
      if digits.len() == 0 {
        return 0;
      }
      panic!("invalid bigint digits len for u32 conversion");
    }
    digits[0]
  }

  pub fn bigint(&self, val: i32) -> BigInt {
    self.modd(&BigInt::from(val))
  }

  pub fn biguint(&self, val: u32) -> BigInt {
    self.modd(&BigInt::new(Sign::Plus, vec!(val)))
  }

  pub fn g(&self) -> BigInt {
    self.g.clone()
  }

  pub fn zero() -> BigInt {
    BigInt::from(0)
  }

  pub fn one() -> BigInt {
    BigInt::from(1)
  }

  pub fn modd(&self, v: &BigInt) -> BigInt {
    if v.sign() == Sign::Minus {
      return (v + &self.p * (1 + (v * -1) / &self.p)) % &self.p;
    }
    v % &self.p
  }

  pub fn add(&self, v1: &BigInt, v2: &BigInt) -> BigInt {
    self.modd(&(v1 + v2))
  }

  pub fn ladd(&self, v1: &BigInt, v2: &BigInt) -> BigInt {
    v1 + v2
  }

  pub fn mul(&self, v1: &BigInt, v2: &BigInt) -> BigInt {
    self.modd(&(v1 * v2))
  }

  pub fn lmul(&self, v1: &BigInt, v2: &BigInt) -> BigInt {
    v1 * v2
  }

  pub fn sub(&self, v1: &BigInt, v2: &BigInt) -> BigInt {
    self.modd(&(v1 - v2))
  }

  pub fn lsub(&self, v1: &BigInt, v2: &BigInt) -> BigInt {
    v1 - v2
  }

  pub fn neg(&self, v: &BigInt) -> BigInt {
    self.mul(v, &(&self.p-&Field::one()))
  }

  pub fn div(&self, v1: &BigInt, v2: &BigInt) -> BigInt {
    self.mul(v1, &self.inv(v2))
  }

  // exponent should always be >= 0
  pub fn exp(&self, v: &BigInt, e: &BigInt) -> BigInt {
    self.modd(v).modpow(e, &self.p)
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
    let mut val = self.modd(&v.clone());
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
    self.modd(&x)
  }

  pub fn random(&self) -> BigInt {
    let mut rng = rand::thread_rng();
    self.bigint(rng.gen())
  }

  pub fn sample(&self, input: &BigInt) -> BigInt {
    self.modd(input)
  }

  pub fn coset(&self, size: u32, offset: &BigInt) -> Vec<BigInt> {
    {
      let cache = self.coset_cache.read().unwrap();
      if let Some(coset) = cache.get(&(size, offset.clone())) {
        return coset.clone();
      }
    }
    let mut cache = self.coset_cache.write().unwrap();
    // double check that it hasn't been built/inserted
    // while we're waiting
    if let Some(coset) = cache.get(&(size, offset.clone())) {
      return coset.clone();
    }
    let generator = self.generator(&BigInt::from(size));
    // build, insert, and return
    let domain = self.domain(&generator, size);
    let mut coset: Vec<BigInt> = Vec::new();
    for v in domain {
      coset.push(self.mul(&v, &offset));
    }
    let d = coset.clone();
    cache.insert((size, offset.clone()), coset);
    d
  }

  // retrieve the full domain, build if necessary
  pub fn domain(&self, generator: &BigInt, size: u32) -> Vec<BigInt> {
    {
      let cache = self.group_cache.read().unwrap();
      if let Some(domain) = cache.get(&(generator.clone(), size)) {
        return domain.clone();
      }
    }
    let mut cache = self.group_cache.write().unwrap();
    // double check that it hasn't been built/inserted
    // while we're waiting
    if let Some(domain) = cache.get(&(generator.clone(), size)) {
      return domain.clone();
    }
    // build, insert, and return
    let mut domain: Vec<BigInt> = Vec::new();
    domain.push(BigInt::from(1));
    for i in 1..size {
      domain.push(self.mul(&domain[domain.len() - 1], &generator));
    }
    let d = domain.clone();
    cache.insert((generator.clone(), size), domain);
    d
  }

}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn should_make_bigint() {
    let p = BigInt::from(101);
    let f = Field::new(p, BigInt::from(0));
    assert_eq!(f.bigint(0), BigInt::new(Sign::Minus, vec!(0)));
    assert_eq!(f.bigint(-29), BigInt::new(Sign::Plus , vec!(72)));
    assert_eq!(f.bigint(32), BigInt::new(Sign::Plus, vec!(32)));
    assert_eq!(f.biguint(0), BigInt::new(Sign::Minus, vec!(0)));
    assert_eq!(f.biguint(32), BigInt::new(Sign::Plus, vec!(32)));
  }

  #[test]
  fn should_add_two_elements() {
    let p = BigInt::from(101);
    let g = BigInt::from(0);
    let f = Field::new(p, g);

    let x = f.bigint(40);
    let y = f.bigint(90);

    assert_eq!(f.add(&x, &y), f.bigint(29));
  }

  #[test]
  fn should_mul_two_elements() {
    let p = BigInt::from(101);
    let g = BigInt::from(0);
    let f = Field::new(p, g);

    let x = f.bigint(40);
    let y = f.bigint(90);

    assert_eq!(f.mul(&x, &y), f.bigint(65));
  }

  #[test]
  fn should_sub_two_elements() {
    let p = BigInt::from(101);
    let g = BigInt::from(0);
    let f = Field::new(p, g);

    let x = f.bigint(2);
    let y = f.bigint(20);

    assert_eq!(f.sub(&x, &y), f.biguint(83));
  }

  #[test]
  fn should_get_generator() {
    let p = BigInt::from(3221225473_u32);
    let f_g = BigInt::from(5);
    let f = Field::new(p, f_g);

    for i in 1..10 {
      let g = f.generator(&f.biguint(u32::pow(2, i)));
      assert_eq!(f.exp(&g, &f.biguint(u32::pow(2, i))), f.bigint(1));
    }
  }

  #[test]
  fn should_get_inverse() {
    let p = BigInt::from(3221225473_u32);
    let f_g = BigInt::from(5);
    let f = Field::new(p, f_g);

    for i in 1..99 {
      let v = f.biguint(i);
      let inv = f.inv(&v);
      assert_eq!(f.mul(&inv, &v), Field::one());
    }
  }
}
