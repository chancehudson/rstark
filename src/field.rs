use num_bigint::{BigInt, Sign};
use num_integer::Integer;
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
  coset_cache: RwLock<HashMap<(u32, BigInt), Vec<BigInt>>>,
  generator_cache: HashMap<u32, (BigInt, BigInt)>
}

impl Field {
  pub fn new(p: BigInt, g: BigInt) -> Field {
    let mut f = Field {
      p,
      g,
      group_cache: RwLock::new(HashMap::new()),
      coset_cache: RwLock::new(HashMap::new()),
      generator_cache: HashMap::new()
    };
    if f.p.bits() <= 32 {
      // we're likely in the 101 field in tests
      // don't build the domain cache
      return f;
    }

    // build a cache of generators and inverted generators
    let mut start = 2;
    let mut sizes: Vec<u32> = Vec::new();
    let mut generators: Vec<BigInt> = Vec::new();
    for _ in 0..31 {
      generators.push(f.generator(&BigInt::from(start)));
      sizes.push(start);
      start *= 2;
    }
    let generators_inv = f.inv_batch(&generators);
    for i in 0..(generators.len()) {
      f.generator_cache.insert(sizes[i], (generators[i].clone(), generators_inv[i].clone()));
    }
    f
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

  pub fn p(&self) -> &BigInt {
    &self.p
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
    self.modd(&(&self.p - v))
  }

  pub fn div(&self, v1: &BigInt, v2: &BigInt) -> BigInt {
    self.mul(v1, &self.inv(v2))
  }

  // exponent should always be >= 0
  pub fn exp(&self, v: &BigInt, e: &BigInt) -> BigInt {
    self.modd(v).modpow(e, &self.p)
  }

  pub fn generator_cache(&self, size: &u32) -> (BigInt, BigInt) {
    if let Some(v) = self.generator_cache.get(size) {
      return v.clone();
    }
    let g = self.generator(&BigInt::from(size.clone()));
    let g_inv = self.inv(&g);
    return (g, g_inv);
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
    let e_gcd = v.clone().extended_gcd(&self.p);
    if e_gcd.gcd != BigInt::from(1) {
      panic!("modinv does not exist");
    }
    self.add(&e_gcd.x, &self.p)
  }

  pub fn inv_batch(&self, values: &Vec<BigInt>) -> Vec<BigInt> {
    let mut last = BigInt::from(1);
    let mut result: Vec<BigInt> = Vec::new();
    for v in values {
      result.push(last.clone());
      last = self.mul(&last, &v);
    }
    last = self.inv(&last);
    for i in (0..values.len()).rev() {
      result[i] = self.mul(&result[i], &last);
      last = self.mul(&values[i], &last);
    }
    result
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
    for _ in 1..size {
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
