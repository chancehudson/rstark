use num_bigint::{BigInt, Sign};
use crate::field::Field;
use std::rc::Rc;

#[derive(Clone)]
pub struct Polynomial {
  field: Rc<Field>,
  coefs: Vec<BigInt>
}

impl Polynomial {
  pub fn field(&self) -> &Rc<Field> {
    &self.field
  }

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
        out[e] = self.field.add(&out[e], &self.field.lmul(&self.coefs[j], &poly.coefs()[i]));
      }
    }
    self.coefs = out;
    self.trim();
    self
  }

  pub fn mul_scalar(& mut self, v: &BigInt) -> &Self {
    for i in 0..self.coefs.len() {
      self.coefs[i] = self.field.mul(v, &self.coefs[i]);
    }
    self
  }

  pub fn exp(& mut self, v: usize) -> &Self {
    let mut out = Polynomial::new(&self.field);
    out.term(&BigInt::from(1), 0);
    for _ in 0..v {
      out.mul(&self);
    }
    self.coefs = out.coefs;
    self
  }

  // if we're scaling the polynomial using a generator point or similar
  // we probably already have a list of the exponents laying around
  pub fn scale_precalc(& mut self, v: &BigInt, exps: &Vec<BigInt>) -> &Self {
    self.coefs = self.coefs.iter().enumerate().map(|(exp, coef)| {
      return self.field.mul(&exps[exp], &coef);
    }).collect();
    self
  }

  pub fn scale(& mut self, v: &BigInt) -> &Self {
    self.coefs = self.coefs.iter().enumerate().map(|(exp, coef)| {
      return self.field.mul(&self.field.exp(&v, &BigInt::from(exp)), &coef);
    }).collect();
    self
  }

  // compose `poly` into `this`
  pub fn compose(& mut self, poly: &Polynomial) -> &Self {
    let mut out = Polynomial::new(&self.field);
    for (exp, coef) in self.coefs.iter().enumerate() {
      let mut p = poly.clone();
      p.exp(exp);
      p.mul_scalar(&coef);
      out.add(&p);
    }
    self.coefs = out.coefs;
    self
  }

  // trim trailing zero coefficient
  pub fn trim(& mut self) {
    let mut new_len = self.coefs.len();
    let zero = Field::zero();
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
  pub fn eval(&self, v: &BigInt) -> BigInt {
    let mut out = self.field.mul(&v, &self.coefs[self.coefs.len() - 1]);
    for coef in self.coefs[1..(self.coefs.len()-1)].iter().rev() {
      out = self.field.mul(&v, &self.field.ladd(&coef, &out));
    }
    out = self.field.add(&out, &self.coefs[0]);
    out
  }

  pub fn eval_batch(&self, vals: &Vec<BigInt>) -> Vec<BigInt> {
    vals.iter().map(|v| self.eval(v)).collect()
  }

  // https://en.wikipedia.org/wiki/Polynomial_evaluation#Multipoint_evaluation
  pub fn eval_batch_fast(&self, vals: &Vec<BigInt>) -> Vec<BigInt> {
    self.eval_batch_fast_slice(&vals[0..])
  }

  pub fn eval_batch_fast_slice(&self, vals: &[BigInt]) -> Vec<BigInt> {
    if vals.len() < 8 {
      return self.eval_batch(&vals.to_vec());
    }
    self.eval_batch_fast_(&vals[0..])
  }

  fn eval_batch_fast_(&self,
    vals: &[BigInt],
  ) -> Vec<BigInt> {
    if vals.len() == 0{
      return Vec::new();
    }
    if vals.len() == 1 {
      return vec!(self.eval(&vals[0]));
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

  pub fn eval_batch_coset(&self, offset: &BigInt, size: u32) -> Vec<BigInt> {
    let mut scaled = self.clone();
    scaled.scale(offset);
    let generator = self.field.generator(&BigInt::from(size));
    let domain = self.field.domain(&generator, size);
    Self::eval_fft(&scaled.coefs(), &domain, &self.field, 1, 0)
  }

  // Evaluate a polynomial over a multiplicative subgroup
  // domain cannot be a coset
  pub fn eval_batch_fft(&self, domain: &Vec<BigInt>) -> Vec<BigInt> {
    Polynomial::eval_fft(&self.coefs(), domain, &self.field, 1, 0)
  }

  pub fn eval_fft(coefs: &Vec<BigInt>, domain: &Vec<BigInt>, field: &Rc<Field>, slice_len: u32, offset: u32) -> Vec<BigInt> {
    let slice_len_usize = usize::try_from(slice_len).unwrap();

    if domain.len()/slice_len_usize == 1 {
      return vec!(coefs.get(usize::try_from(offset).unwrap()).unwrap_or(&BigInt::from(0)).clone());
    }
    let left_out = Self::eval_fft(coefs, domain, field, slice_len*2, offset);
    let right_out = Self::eval_fft(coefs, domain, field, slice_len*2, offset + slice_len);

    let mut out1: Vec<BigInt> = Vec::new();
    let mut out2: Vec<BigInt> = Vec::new();

    for i in 0..(left_out.len()) {
      let x = &left_out[i];
      let y = &right_out[i];
      let y_root = field.mul(y, &domain[i*slice_len_usize]);
      // bring the values into the field using simple arithmetic
      // instead of relying on modulus
      // offers a small speedup
      let mut o1 = x + &y_root;
      if &o1 > field.p() {
        o1 -= field.p();
      }
      out1.push(o1);
      let mut o2 = x - &y_root;
      if o2 < Field::zero() {
        o2 += field.p();
      }
      out2.push(o2);
    }
    out1.extend(out2);
    out1
  }

  pub fn eval_fft_inv(coefs: &Vec<BigInt>, domain_inv: &Vec<BigInt>, field: &Rc<Field>) -> Vec<BigInt> {
    if coefs.len() == 1 {
      return vec!(coefs[0].clone());
    }
    let len_inv = field.inv(&BigInt::from(u32::try_from(coefs.len()).unwrap()));
    let out = Self::eval_fft(coefs, &domain_inv, field, 1, 0);
    out.iter().map(|v| field.mul(&v, &len_inv)).collect()
  }

  pub fn mul_fft(poly1: &Polynomial, poly2: &Polynomial, field: &Rc<Field>) -> Polynomial {
    if poly1.degree() + poly2.degree() < 20 {
      let mut o = poly1.clone();
      o.mul(&poly2);
      return o;
    }
    let out_degree = 2*std::cmp::max(poly1.degree(), poly2.degree());
    let domain_size = 2_u32.pow(u32::try_from(out_degree).unwrap().ilog2() + 1);
    let generator = field.generator(&BigInt::from(domain_size));
    let domain = field.domain(&generator, domain_size);

    let x1 = Self::eval_fft(poly1.coefs(), &domain, field, 1, 0);
    let x2 = Self::eval_fft(poly2.coefs(), &domain, field, 1, 0);

    let mut x3: Vec<BigInt> = Vec::new();
    for i in 0..domain.len() {
      x3.push(field.mul(&x1[i], &x2[i]));
    }
    let generator_inv = field.inv(&generator);
    let domain_inv = field.domain(&generator_inv, domain_size);

    let out = Self::eval_fft_inv(&x3, &domain_inv, field);
    let mut p = Polynomial {
      field: Rc::clone(field),
      coefs: out
    };
    p.trim();
    p
  }

  pub fn div_coset(
    poly1: &Polynomial,
    poly2: &Polynomial,
    offset: &BigInt,
    generator: &BigInt,
    size: u32,
    field: &Rc<Field>
  ) -> Polynomial {
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
    let offset_domain = field.domain(&offset, order);

    let mut poly1_scaled = poly1.clone();
    poly1_scaled.scale_precalc(&offset, &offset_domain);
    let mut poly2_scaled = poly2.clone();
    poly2_scaled.scale_precalc(&offset, &offset_domain);

    let poly1_codeword = Self::eval_fft(poly1_scaled.coefs(), &domain, field, 1, 0);
    let poly2_codeword = Self::eval_fft(poly2_scaled.coefs(), &domain, field, 1, 0);

    let out = poly1_codeword.iter().enumerate().map(|(i, val)| {
      field.div(&val, &poly2_codeword[i])
    }).collect();

    let scaled_coefs = Self::eval_fft_inv(&out, &domain_inv, field);
    let mut scaled_poly = Polynomial::new(field);
    for i in 0..(poly1.degree() - poly2.degree() + 1) {
      scaled_poly.term(&scaled_coefs[i], u32::try_from(i).unwrap());
    }
    scaled_poly.scale(&field.inv(&offset));
    scaled_poly
  }

  // remove and return the largest non-zero coefficient
  // coef, exp
  pub fn pop_term(& mut self) -> (BigInt, usize) {
    let zero = Field::zero();
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

  pub fn safe_div(&self, divisor: &Polynomial) -> Polynomial {
    let (q, r) = self.div(divisor);
    if !r.is_zero() {
      panic!("non-zero remainder in division");
    }
    q
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
      // inter.sub(&Polynomial::mul_fft(&t, &divisor, &self.field));
    }
    (q, inter)
  }

  pub fn lagrange(x_vals: &Vec<BigInt>, y_vals: &Vec<BigInt>, field: &Rc<Field>) -> Polynomial {
    if x_vals.len() != y_vals.len() {
      panic!("lagrange mismatch x/y array length");
    }
    let mut numerator = Polynomial::new(field);
    numerator.term(&field.bigint(1), 0);
    for v in x_vals {
      let mut poly = Polynomial::new(field);
      poly.term(&field.bigint(1), 1);
      poly.term(&field.neg(&v), 0);
      numerator.mul(&poly);
    }
    let mut polynomials: Vec<Polynomial> = Vec::new();
    for i in 0..x_vals.len() {
      let mut denominator = Field::one();
      for j in 0..x_vals.len() {
        if i == j {
          continue;
        }
        denominator = field.mul(&denominator, &(&x_vals[i] - &x_vals[j]));
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

  pub fn interpolate_fft(x_vals: &Vec<BigInt>, y_vals: &Vec<BigInt>, field: &Rc<Field>) -> Polynomial {
    Self::interpolate_fft_slice(&x_vals[0..], &y_vals[0..], field)
  }

  pub fn interpolate_fft_slice(x_vals: &[BigInt], y_vals: &[BigInt], field: &Rc<Field>) -> Polynomial {
    if x_vals.len() != y_vals.len() {
      panic!("x/y len mismatch");
    }
    if x_vals.len() == 0 {
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

    let left_targets = y_vals[0..half].iter().enumerate().map(|(i, v)| field.div(v, &left_offset[i])).collect::<Vec<BigInt>>();
    let right_targets = y_vals[half..].iter().enumerate().map(|(i, v)| field.div(v, &right_offset[i])).collect::<Vec<BigInt>>();

    let mut left_interpolant = Self::interpolate_fft_slice(&x_vals[0..half], &left_targets[0..], field);
    let mut right_interpolant = Self::interpolate_fft_slice(&x_vals[half..], &right_targets[0..], field);

    left_interpolant.mul(&right_zeroifier);
    right_interpolant.mul(&left_zeroifier);
    left_interpolant.add(&right_interpolant);
    left_interpolant
  }

  pub fn interpolate_fft_batch(x_vals: &[BigInt], y_vals: &[Vec<BigInt>], field: &Rc<Field>) -> Vec<Polynomial> {
    if x_vals.len() == 0 {
      return vec!(Polynomial::new(field));
    }
    if x_vals.len() == 1 {
      return y_vals.iter().map(|vals| {
        let mut p = Polynomial::new(field);
        p.term(&vals[0], 0);
        p
      }).collect();
    }
    let half = x_vals.len() >> 1;

    let left_zeroifier = Self::zeroifier_fft_slice(&x_vals[0..half], field);
    let right_zeroifier = Self::zeroifier_fft_slice(&x_vals[half..], field);

    let left_offset = right_zeroifier.eval_batch_fast_slice(&x_vals[0..half]);
    let right_offset = left_zeroifier.eval_batch_fast_slice(&x_vals[half..]);

    let left_targets: Vec<Vec<BigInt>> = y_vals.iter().map(|vals| {
      vals[0..half].iter().enumerate().map(|(i, v)| field.div(v, &left_offset[i])).collect::<Vec<BigInt>>()
    }).collect();
    let right_targets: Vec<Vec<BigInt>> = y_vals.iter().map(|vals| {
      vals[half..].iter().enumerate().map(|(i, v)| field.div(v, &right_offset[i])).collect::<Vec<BigInt>>()
    }).collect();

    let left_interpolant = Self::interpolate_fft_batch(&x_vals[0..half], &left_targets[0..], field);
    let mut right_interpolant = Self::interpolate_fft_batch(&x_vals[half..], &right_targets[0..], field);

    left_interpolant.iter().enumerate().map(|(i, poly)| {
      let mut left = poly.clone();
      left.mul(&right_zeroifier);
      right_interpolant[i].mul(&left_zeroifier);
      left.add(&right_interpolant[i]);
      left
    }).collect()
  }

  pub fn test_colinearity(x_vals: &Vec<BigInt>, y_vals: &Vec<BigInt>, field: &Rc<Field>) -> bool {
    let poly = Polynomial::lagrange(x_vals, y_vals, field);
    return poly.degree() <= 1;
  }

  pub fn zeroifier_fft(points: &Vec<BigInt>, field: &Rc<Field>) -> Polynomial {
    Self::zeroifier_fft_slice(&points[0..], field)
  }

  pub fn zeroifier_fft_slice(points: &[BigInt], field: &Rc<Field>) -> Polynomial {
    if points.len() == 0 {
      return Polynomial::new(field);
    }
    if points.len() == 1 {
      let mut p = Polynomial::new(field);
      p.term(&field.neg(&points[0]), 0);
      p.term(&BigInt::from(1), 1);
      return p;
    }
    let half = points.len() >> 1;
    let left = Self::zeroifier_fft_slice(&points[0..(half)], field);
    let right = Self::zeroifier_fft_slice(&points[(half)..], field);
    return Self::mul_fft(&left, &right, field);
  }

  pub fn zeroifier(points: &Vec<BigInt>, field: &Rc<Field>) -> Polynomial {
    Self::zeroifier_slice(&points[0..], field)
  }

  pub fn zeroifier_slice(points: &[BigInt], field: &Rc<Field>) -> Polynomial {
    let mut out = Polynomial::new(field);
    out.term(&BigInt::from(1), 0);
    let mut x = Polynomial::new(field);
    x.term(&BigInt::from(1), 1);
    for p in points {
      out.mul(&x.clone().term(&field.neg(&p), 0));
    }
    out
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn should_test_colinearity() {
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    poly.term(&BigInt::from(-9), 0);
    poly.term(&BigInt::from(2), 1);
    let mut x_vals: Vec<BigInt> = Vec::new();
    let mut y_vals: Vec<BigInt> = Vec::new();
    for i in 0..3 {
      x_vals.push(BigInt::from(i));
      y_vals.push(poly.eval(&BigInt::from(i)));
    }
    assert!(Polynomial::test_colinearity(&x_vals, &y_vals, &f));
  }

  #[test]
  fn should_compose_polynomial() {
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
    let f = Rc::new(Field::new(p, g));

    let mut root = Polynomial::new(&f);
    root.term(&BigInt::from(99), 0);
    root.term(&BigInt::from(2), 1);
    root.term(&BigInt::from(4), 2);

    let mut inpoly = Polynomial::new(&f);
    inpoly.term(&BigInt::from(2), 2);
    inpoly.term(&BigInt::from(12), 0);

    let mut expected = Polynomial::new(&f);
    expected.term(&BigInt::from(99), 0);
    expected.add(&inpoly.clone().mul_scalar(&BigInt::from(2)));
    {
      let mut i = inpoly.clone();
      i.exp(2);
      i.mul_scalar(&BigInt::from(4));
      expected.add(&i);
    }

    assert!(root.compose(&inpoly).is_equal(&expected));
  }

  #[test]
  fn should_exp_polynomial() {
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    poly.term(&BigInt::from(2), 0);

    for i in 0..10 {
      let mut expected = Polynomial::new(&f);
      expected.term(&f.exp(&BigInt::from(2), &BigInt::from(i)), 0);
      assert!(poly.clone().exp(i).is_equal(&expected));
    }
  }

  #[test]
  fn should_scale_polynomial() {
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
    let f = Rc::new(Field::new(p, g));

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

    assert!(expected.is_equal(poly.scale(&BigInt::from(2))));
  }

  #[test]
  fn should_interpolate_lagrange() {
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
    let f = Rc::new(Field::new(p, g));

    let size = 32;
    let gen = f.generator(&f.bigint(size));
    let mut x_vals: Vec<BigInt> = Vec::new();
    let mut y_vals: Vec<BigInt> = Vec::new();
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
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    poly.term(&f.bigint(1), 0);

    let mut divisor = Polynomial::new(&f);
    divisor.term(&f.bigint(0), 3);
    poly.div(&divisor);
  }

  #[test]
  fn should_divide_polynomial() {
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
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
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
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
  fn should_eval_polynomial_with_batch_fast() {
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    for i in 0..50 {
      poly.term(&f.random(), i);
    }

    let size = 2_u32.pow(7);
    let mut G: Vec<BigInt> = Vec::new();
    for i in 0..size {
      G.push(BigInt::from(i));
    }

    let actual = poly.eval_batch(&G);
    let out = poly.eval_batch_fast(&G);
    for i in 0..usize::try_from(size).unwrap() {
      assert_eq!(actual[i], out[i]);
    }
  }

  #[test]
  fn should_eval_polynomial_with_fft() {
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    for i in 0..50 {
      poly.term(&f.random(), i);
    }

    let size = 2_u32.pow(8);
    let sub_g = f.generator(&BigInt::from(size));
    let mut G: Vec<BigInt> = Vec::new();
    for i in 0..size {
      G.push(f.exp(&sub_g, &BigInt::from(i)));
    }

    let actual = poly.eval_batch(&G);
    let out = poly.eval_batch_fft(&G);
    for i in 0..usize::try_from(size).unwrap() {
      assert_eq!(actual[i], out[i]);
    }
  }

  #[test]
  fn should_check_polynomial_equality() {
    let p = BigInt::from(101);
    let g = BigInt::from(0);
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
    let p = BigInt::from(101);
    let g = BigInt::from(0);
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
    let p = BigInt::from(101);
    let g = BigInt::from(0);
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
    let p = BigInt::from(101);
    let g = BigInt::from(0);
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

  #[test]
  fn should_multiply_polynomials_fft() {
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
    let f = Rc::new(Field::new(p, g));

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
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
    let f = Rc::new(Field::new(p, g));

    let s = 128;
    let size = BigInt::from(s);
    let domain_g = f.generator(&size);
    let mut domain: Vec<BigInt> = Vec::new();
    for i in 0..s {
      domain.push(f.exp(&domain_g, &BigInt::from(i)));
    }

    let poly = Polynomial::zeroifier(&domain, &f);
    for i in 0..s {
      assert_eq!(poly.eval(&domain[i]), BigInt::from(0));
    }

    let mut is_zero = true;
    for i in poly.coefs() {
      if i != &BigInt::from(0) {
        is_zero = false;
      }
    }
    assert!(!is_zero);
  }

  #[test]
  fn should_build_zeroifier_polynomial_domain() {
    let p = BigInt::from(3221225473_u32);
    let g = BigInt::from(5);
    let f = Rc::new(Field::new(p, g));

    let size = 128_u32;
    let generator = f.generator(&BigInt::from(size));
    let domain = f.domain(&generator, size);
    let zeroifier = Polynomial::zeroifier(&domain, &f);
    let mut zeroifier_fft = Polynomial::zeroifier_fft(&domain, &f);
    zeroifier_fft.trim();

    for v in domain {
      assert_eq!(BigInt::from(0), zeroifier_fft.eval(&v));
      assert_eq!(BigInt::from(0), zeroifier.eval(&v));
    }

    assert!(zeroifier.is_equal(&zeroifier_fft));
  }
}
