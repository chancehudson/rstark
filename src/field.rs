#[cfg(test)]
use num_bigint::BigInt;

use rand::Rng;
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::sync::RwLock;

static upper: u128 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF << 64_u128;

#[derive(Serialize, Deserialize)]
pub struct Field {
  p: u128,
  g: u128,
  // generator, size size keyed to contents
  group_cache: RwLock<HashMap<(u128, u32), Vec<u128>>>,
  // group size, offset keyed to contents
  coset_cache: RwLock<HashMap<(u32, u128), Vec<u128>>>,
  generator_cache: HashMap<u32, (u128, u128)>
}

impl Field {
  pub fn new(p: u128, g: u128) -> Field {
    let mut f = Field {
      p,
      g,
      group_cache: RwLock::new(HashMap::new()),
      coset_cache: RwLock::new(HashMap::new()),
      generator_cache: HashMap::new()
    };
    if p.ilog2() < 32 {
      // we're likely in the 101 field in tests
      // don't build the domain cache
      return f;
    }

    // build a cache of generators and inverted generators
    let mut start: u32 = 1;
    let mut sizes: Vec<u32> = Vec::new();
    let mut generators: Vec<u128> = Vec::new();
    for _ in 0..31 {
      start *= 2;
      generators.push(f.generator(&start));
      sizes.push(start);
    }
    let generators_inv = f.inv_batch(&generators);
    for i in 0..(generators.len()) {
      f.generator_cache.insert(sizes[i], (generators[i].clone(), generators_inv[i].clone()));
    }
    f
  }

  pub fn bigint_to_u32(v: &u128) -> u32 {
    if v.ilog2() >= 32 {
      panic!("invalid bigint digits len for u32 conversion");
    }
    u32::try_from(v.clone()).unwrap()
  }

  pub fn p(&self) -> &u128 {
    &self.p
  }

  pub fn g(&self) -> u128 {
    self.g.clone()
  }

  pub fn zero() -> u128 {
    0
  }

  pub fn one() -> u128 {
    1
  }

  pub fn modd(&self, v: &u128) -> u128 {
    v % &self.p
  }

  pub fn add(&self, v1: &u128, v2: &u128) -> u128 {
    let (r, overflowed) = v1.overflowing_add(v2.clone());
    if !overflowed {
      return r % self.p;
    }
    return r + (u128::MAX - self.p) + 1;
  }


  // multiply two numbers modulo self.p
  pub fn mul(&self, v1: &u128, v2: &u128) -> u128 {
    if ((v1|v2) & upper) == 0 {
      return (v1 * v2) % self.p;
    }

    let mut a;
    let mut b;
    if v2 < v1 {
      a = v2.clone();
      b = v1.clone();
    } else {
      a = v1.clone();
      b = v2.clone();
    }
    let mut r = 0_u128;

    while a != 0 {
      if a & 1 == 1 {
        if b >= self.p - r {
          r -= self.p - b;
        } else {
          r += b;
        }
      }
      a >>= 1;
      if b >= self.p - b {
        b -= self.p - b;
      } else {
        b += b;
      }
    }

    r
}

  pub fn sub(&self, v1: &u128, v2: &u128) -> u128 {
    if v1 >= v2 {
      return v1 - v2;
    }
    return v1 + (self.p - v2);
  }

  pub fn neg(&self, v: &u128) -> u128 {
    &self.p - v
  }

  pub fn div(&self, v1: &u128, v2: &u128) -> u128 {
    self.mul(v1, &self.inv(v2))
  }

  // calculate v^e using 2^k-ary method
  // https://en.wikipedia.org/wiki/Exponentiation_by_squaring
  pub fn exp(&self, v: &u128, e: &u128) -> u128 {
    if e == &0 {
      return 1;
    }
    if e == &1 {
      return v.clone();
    }
    if e == &2 {
      return self.mul(v, v);
    }
    let mut e_ = e.clone();
    let mut v_ = v.clone();
    let mut t = 1_u128;
    while e_ > 0 {
      if e_ % 2 != 0 {
        t = self.mul(&t, &v_);
      }
      v_ = self.mul(&v_, &v_);
      e_ = e_ >> 1;
    }
    t % self.p
  }

  pub fn generator_cache(&self, size: &u32) -> (u128, u128) {
    if let Some(v) = self.generator_cache.get(size) {
      return v.clone();
    }
    let g = self.generator(size);
    let g_inv = self.inv(&g);
    return (g, g_inv);
  }

  pub fn generator(&self, size: &u32) -> u128 {
    if size.ilog2() >= self.p.ilog2() {
      panic!("requested subgroup is larger than field");
    }
    let numer = &self.p - Field::one();
    let size_128 = u128::try_from(size.clone()).unwrap();
    let exp = &numer / size_128;
    if exp.clone() * size_128 != numer {
      panic!("subgroup is not a divisor of field");
    }
    self.exp(&self.g, &exp)
  }

  pub fn inv(&self, v: &u128) -> u128 {
    let mut val = self.modd(v);
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
      y = self.sub(&x, &self.mul(&q, &y));
      x = t;
    }
    self.modd(&x)
  }

  pub fn inv_batch(&self, values: &Vec<u128>) -> Vec<u128> {
    let mut last = 1_u128;
    let mut result: Vec<u128> = Vec::new();
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

  pub fn random(&self) -> u128 {
    let mut rng = rand::thread_rng();
    self.modd(&rng.gen::<u128>())
  }

  pub fn sample(&self, input: &u128) -> u128 {
    self.modd(input)
  }

  pub fn sample_bytes(&self, input: &[u8]) -> u128 {
    let mut out = 0_u128;
    for i in 0..16 {
      let i_b = u128::try_from(i).unwrap();
      out = self.add(&out, &self.mul(&u128::try_from(input[i]).unwrap(), &self.exp(&2, &i_b)));
    }
    out
  }

  pub fn coset(&self, size: u32, offset: &u128) -> Vec<u128> {
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
    let generator = self.generator(&size);
    // build, insert, and return
    let domain = self.domain(&generator, size);
    let mut coset: Vec<u128> = Vec::new();
    for v in domain {
      coset.push(self.mul(&v, &offset));
    }
    let d = coset.clone();
    cache.insert((size, offset.clone()), coset);
    d
  }

  // retrieve the full domain, build if necessary
  pub fn domain(&self, generator: &u128, size: u32) -> Vec<u128> {
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
    let mut domain: Vec<u128> = Vec::new();
    domain.push(Field::one());
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
    let p = 101_u128;
    let f = Field::new(p, 0_u128);
    assert_eq!(0, 0);
    assert_eq!(f.neg(&29), 72);
    assert_eq!(32, 32);
  }

  #[test]
  fn should_add_two_elements() {
    let p = 101_u128;
    let g = Field::zero();
    let f = Field::new(p, g);

    let x = 40;
    let y = 90;

    assert_eq!(f.add(&x, &y), 29);
  }

  #[test]
  fn should_mul_two_elements() {
    {
      let p = 101_u128;
      let g = Field::zero();
      let f = Field::new(p, g);

      let x = 40;
      let y = 90;

      assert_eq!(f.mul(&x, &y), 65);
    }
    {
      let p = 3221225473_u128;
      let f_g = 5_u128;
      let f = Field::new(p, f_g);

      for _ in 0..100 {
        let a = f.random();
        let b = f.random();
        let expected = (a*b) % p;
        assert_eq!(f.mul(&a, &b), expected);
      }
    }
    {
      let p = 1 + 407 * 2_u128.pow(119);
      let g = 85408008396924667383611388730472331217_u128;
      let f = Field::new(p, g.clone());
      for _ in 0..1000 {
        let a = f.random();
        let b = f.random();
        let expected = BigInt::from(a)*BigInt::from(b) % BigInt::from(p);
        assert_eq!(BigInt::from(f.mul(&a, &b)), expected);
      }
    }
  }

  #[test]
  fn should_sub_two_elements() {
    let p = 101_u128;
    let g = Field::zero();
    let f = Field::new(p, g);

    let x = 2;
    let y = 20;

    assert_eq!(f.sub(&x, &y), 83);
  }

  #[test]
  fn should_get_generator() {
    let p = 3221225473_u128;
    let f_g = 5_u128;
    let f = Field::new(p, f_g);

    for i in 1..10 {
      let g = f.generator(&u32::pow(2, i));
      assert_eq!(f.exp(&g, &u128::pow(2, i)), 1);
    }
  }

  #[test]
  fn should_get_inverse() {
    let p = 3221225473_u128;
    let f_g = 5_u128;
    let f = Field::new(p, f_g);

    for i in 1..99 {
      let v = f.random();
      let inv = f.inv(&v);
      assert_eq!(f.mul(&inv, &v), Field::one());
    }
  }

  #[test]
  fn should_exponentiate() {
    let p = 3221225473_u128;
    let f_g = 5_u128;
    let f = Field::new(p, f_g);

    for _ in 0..100 {
      let v = f.random();
      let mut exp = 1_u128;
      for i in 0..50_u128 {
        assert_eq!(exp, f.exp(&v, &i));
        exp = f.mul(&exp, &v);
      }
    }
  }
}
