
use crate::field::Field;
use std::rc::Rc;

#[derive(Clone)]
pub struct Polynomial {
  field: Rc<Field>,
  coefs: Vec<u128>
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

  pub fn coefs(&self) -> &Vec<u128> {
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

  pub fn term(& mut self, coef: &u128, exp: u32) -> &Self {
    let s = usize::try_from(exp).unwrap();
    // expand to one longer to handle 0 exponents
    if self.coefs.len() < s+1 {
      self.coefs.resize(s+1, Field::zero());
    }
    self.coefs[s] = self.field.add(&self.coefs[s], &self.field.modd(coef));
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

  // a fast method for multiplying by a single term polynomial
  // with a coefficient of 1
  // e.g. multiplying by x^5
  pub fn shift_and_clone(&self, degree: u32) -> Self {
    if degree == 0 {
      return self.clone();
    }
    let degree_usize = usize::try_from(degree).unwrap();
    let mut shifted_coefs = vec!(0_u128; degree_usize);
    shifted_coefs.extend(self.coefs.clone());
    let mut out = Polynomial::new(&self.field);
    out.coefs = shifted_coefs;
    out
  }

  pub fn mul(& mut self, poly: &Polynomial) -> &Self {
    let mut out: Vec<u128> = Vec::new();
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

  pub fn mul_scalar(& mut self, v: &u128) -> &Self {
    if v == &0 {
      self.coefs = vec!();
      return self;
    }
    for i in 0..self.coefs.len() {
      self.coefs[i] = self.field.mul(v, &self.coefs[i]);
    }
    self
  }

  pub fn exp(& mut self, v: usize) -> &Self {
    let mut out = Polynomial::new(&self.field);
    out.term(&1, 0);
    for _ in 0..v {
      out.mul(&self);
    }
    self.coefs = out.coefs;
    self
  }

  // if we're scaling the polynomial using a generator point or similar
  // we probably already have a list of the exponents laying around
  pub fn scale_precalc(& mut self, _v: &u128, exps: &Vec<u128>) -> &Self {
    self.coefs = self.coefs.iter().enumerate().map(|(exp, coef)| {
      return self.field.mul(&exps[exp], &coef);
    }).collect();
    self
  }

  pub fn scale(& mut self, v: &u128) -> &Self {
    self.coefs = self.coefs.iter().enumerate().map(|(exp, coef)| {
      return self.field.mul(&self.field.exp(&v, &u128::try_from(exp).unwrap()), &coef);
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
  pub fn eval(&self, v: &u128) -> u128 {
    if self.coefs.len() == 0 {
      return 0;
    }
    if self.coefs.len() == 1 {
      return self.coefs[0].clone();
    }
    let mut out = self.field.mul(&v, &self.coefs[self.coefs.len() - 1]);
    for coef in self.coefs[1..(self.coefs.len()-1)].iter().rev() {
      out = self.field.mul(&v, &self.field.add(&coef, &out));
    }
    out = self.field.add(&out, &self.coefs[0]);
    out
  }

  pub fn eval_batch(&self, vals: &Vec<u128>) -> Vec<u128> {
    vals.iter().map(|v| self.eval(v)).collect()
  }

  // https://en.wikipedia.org/wiki/Polynomial_evaluation#Multipoint_evaluation
  pub fn eval_batch_fast(&self, vals: &Vec<u128>) -> Vec<u128> {
    self.eval_batch_fast_slice(&vals[0..])
  }

  pub fn eval_batch_fast_slice(&self, vals: &[u128]) -> Vec<u128> {
    if vals.len() < 8 {
      return self.eval_batch(&vals.to_vec());
    }
    self.eval_batch_fast_(&vals[0..])
  }

  fn eval_batch_fast_(&self,
    vals: &[u128],
  ) -> Vec<u128> {
    if vals.len() == 0 {
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

  pub fn eval_batch_coset(&self, offset: &u128, size: u32) -> Vec<u128> {
    let mut scaled = self.clone();
    let offset_domain = self.field.domain(&offset, size);
    scaled.scale_precalc(offset, &offset_domain);
    let (generator, _) = self.field.generator_cache(&size);
    let domain = self.field.domain(&generator, size);
    Self::eval_fft(&scaled.coefs(), &domain, &self.field)
  }

  // Evaluate a polynomial over a multiplicative subgroup
  // domain cannot be a coset
  pub fn eval_batch_fft(&self, domain: &Vec<u128>) -> Vec<u128> {
    Polynomial::eval_fft(&self.coefs(), domain, &self.field)
  }

  pub fn eval_fft(coefs: &Vec<u128>, domain: &Vec<u128>, field: &Rc<Field>) -> Vec<u128> {
    let mut out = vec!(0_u128; domain.len());
    Self::eval_fft_(coefs, domain, field, 1, 0, 0, domain.len()/2, & mut out);
    out
  }

  // fft implemented by mutating a single `out` vec
  // this saves the cost of alloc/freeing vectors
  // during recursion
  //
  // measured improvement of ~10% compared to previous
  // approach
  pub fn eval_fft_(
    coefs: &Vec<u128>,
    domain: &Vec<u128>,
    field: &Rc<Field>,
    slice_len: usize,
    offset: usize,
    left_dest: usize,
    right_dest: usize,
    out: & mut Vec<u128>
  ) {
    let out_size = domain.len() / slice_len;

    if out_size == 1 {
      if let Some(v) = coefs.get(offset) {
        out[left_dest] = v.clone();
      }
      // otherwise the coef is 0 which is the default value in `out`
      return;
    }

    Self::eval_fft_(coefs, domain, field, slice_len*2, offset, left_dest, left_dest + out_size/4, out);
    Self::eval_fft_(coefs, domain, field, slice_len*2, offset + slice_len, right_dest, right_dest + out_size/4, out);

    for i in 0..(out_size/2) {
      let left_out;
      let right_out;
      {
        let x = &(&*out)[left_dest + i];
        let y = &(&*out)[right_dest + i];
        let y_root = field.mul(y, &domain[i*slice_len]);
        // bring the values into the field using simple arithmetic
        // instead of relying on modulus
        // offers a small speedup
        left_out = field.add(&x, &y_root);
        right_out = field.sub(&x, &y_root);
      }
      out[left_dest + i] = left_out;
      out[right_dest + i] = right_out;
    }
  }

  pub fn eval_fft_inv(coefs: &Vec<u128>, domain_inv: &Vec<u128>, field: &Rc<Field>) -> Vec<u128> {
    if coefs.len() == 1 {
      return vec!(coefs[0].clone());
    }
    let len_inv = field.inv(&u128::try_from(coefs.len()).unwrap());
    let out = Self::eval_fft(coefs, &domain_inv, field);
    out.iter().map(|v| field.mul(&v, &len_inv)).collect()
  }

  pub fn mul_fft(poly1: &Polynomial, poly2: &Polynomial, field: &Rc<Field>) -> Polynomial {
    if poly1.degree() + poly2.degree() < 16 {
      let mut o = poly1.clone();
      o.mul(&poly2);
      return o;
    }
    let out_degree = 2*std::cmp::max(poly1.degree(), poly2.degree());
    let domain_size = 2_u32.pow(u32::try_from(out_degree).unwrap().ilog2() + 1);
    let (generator, generator_inv) = field.generator_cache(&domain_size);
    let domain = field.domain(&generator, domain_size);

    let x1 = Self::eval_fft(poly1.coefs(), &domain, field);
    let x2 = Self::eval_fft(poly2.coefs(), &domain, field);

    let mut x3: Vec<u128> = Vec::new();
    for i in 0..domain.len() {
      x3.push(field.mul(&x1[i], &x2[i]));
    }
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
    offset: &u128,
    generator: &u128,
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

    let poly1_codeword = Self::eval_fft(poly1_scaled.coefs(), &domain, field);
    let poly2_codeword = Self::eval_fft(poly2_scaled.coefs(), &domain, field);

    let poly2_codeword_inv = field.inv_batch(&poly2_codeword);

    let out = poly1_codeword.iter().enumerate().map(|(i, val)| {
      field.mul(&val, &poly2_codeword_inv[i])
    }).collect();

    let scaled_coefs = Self::eval_fft_inv(&out, &domain_inv, field);
    let mut scaled_poly = Polynomial::new(field);
    for i in 0..(poly1.degree() - poly2.degree() + 1) {
      scaled_poly.term(&scaled_coefs[i], u32::try_from(i).unwrap());
    }
    let offset_inv = field.inv(&offset);
    let offset_domain_inv = field.domain(&offset_inv, order);
    scaled_poly.scale_precalc(&offset_inv, &offset_domain_inv);
    scaled_poly
  }

  // remove and return the largest non-zero coefficient
  // coef, exp
  pub fn pop_term(& mut self) -> (u128, usize) {
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
      let mut t = divisor.shift_and_clone(u32::try_from(new_exp).unwrap());
      t.mul_scalar(&new_coef);
      inter.sub(&t);
    }
    inter.trim();
    (q, inter)
  }

  pub fn lagrange(x_vals: &Vec<u128>, y_vals: &Vec<u128>, field: &Rc<Field>) -> Polynomial {
    if x_vals.len() != y_vals.len() {
      panic!("lagrange mismatch x/y array length");
    }
    let mut numerator = Polynomial::new(field);
    numerator.term(&1, 0);
    for v in x_vals {
      let mut poly = Polynomial::new(field);
      poly.term(&1, 1);
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
        denominator = field.mul(&denominator, &field.sub(&x_vals[i], &x_vals[j]));
      }
      let mut n = Polynomial::new(field);
      n.term(&1, 1);
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

  pub fn interpolate_fft(x_vals: &Vec<u128>, y_vals: &Vec<u128>, field: &Rc<Field>) -> Polynomial {
    Self::interpolate_fft_slice(&x_vals[0..], &y_vals[0..], field)
  }

  pub fn interpolate_fft_slice(x_vals: &[u128], y_vals: &[u128], field: &Rc<Field>) -> Polynomial {
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

    let left_offset_invs = field.inv_batch(&left_offset);
    let right_offset_invs = field.inv_batch(&right_offset);

    let left_targets = y_vals[0..half].iter().enumerate().map(|(i, v)| field.mul(v, &left_offset_invs[i])).collect::<Vec<u128>>();
    let right_targets = y_vals[half..].iter().enumerate().map(|(i, v)| field.mul(v, &right_offset_invs[i])).collect::<Vec<u128>>();

    let mut left_interpolant = Self::interpolate_fft_slice(&x_vals[0..half], &left_targets[0..], field);
    let mut right_interpolant = Self::interpolate_fft_slice(&x_vals[half..], &right_targets[0..], field);

    left_interpolant.mul(&right_zeroifier);
    right_interpolant.mul(&left_zeroifier);
    left_interpolant.add(&right_interpolant);
    left_interpolant
  }

  pub fn interpolate_fft_batch(x_vals: &[u128], y_vals: &[Vec<u128>], field: &Rc<Field>) -> Vec<Polynomial> {
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

    let left_offset_invs: Vec<u128> = left_offset.iter().map(|v| field.inv(&v)).collect();
    let right_offset_invs: Vec<u128> = right_offset.iter().map(|v| field.inv(&v)).collect();

    let left_targets: Vec<Vec<u128>> = y_vals.iter().map(|vals| {
      vals[0..half].iter().enumerate().map(|(i, v)| field.mul(v, &left_offset_invs[i])).collect::<Vec<u128>>()
    }).collect();
    let right_targets: Vec<Vec<u128>> = y_vals.iter().map(|vals| {
      vals[half..].iter().enumerate().map(|(i, v)| field.mul(v, &right_offset_invs[i])).collect::<Vec<u128>>()
    }).collect();

    let left_interpolant = Self::interpolate_fft_batch(&x_vals[0..half], &left_targets[0..], field);
    let right_interpolant = Self::interpolate_fft_batch(&x_vals[half..], &right_targets[0..], field);

    let mut out: Vec<Polynomial> = Vec::new();
    for i in 0..left_interpolant.len() {
      let mut left = Polynomial::mul_fft(&left_interpolant[i], &right_zeroifier, &field);
      let right = Polynomial::mul_fft(&right_interpolant[i], &left_zeroifier, &field);
      left.add(&right);
      out.push(left);
    }
    out
  }

  pub fn test_colinearity_batch(x_vals_arr: &Vec<Vec<u128>>, y_vals_arr: &Vec<Vec<u128>>, field: &Rc<Field>) -> bool {
    let mut to_inv: Vec<u128> = Vec::new();
    for x_vals in x_vals_arr {
      to_inv.push(field.sub(&x_vals[1], &x_vals[0]));
      to_inv.push(field.sub(&x_vals[2], &x_vals[1]));
    }
    let inverted = field.inv_batch(&to_inv);
    for (i, y_vals) in y_vals_arr.iter().enumerate() {
      let x_diff_inv_1 = &inverted[2*i];
      let x_diff_inv_2 = &inverted[2*i + 1];
      let y_diff_1 = field.sub(&y_vals[1], &y_vals[0]);
      let y_diff_2 = field.sub(&y_vals[2], &y_vals[1]);
      let slope_1 = field.mul(&y_diff_1, x_diff_inv_1);
      let slope_2 = field.mul(&y_diff_2, x_diff_inv_2);
      if slope_1 != slope_2 {
        return false;
      }
    }
    return true;
  }

  pub fn test_colinearity(x_vals: &Vec<u128>, y_vals: &Vec<u128>, field: &Rc<Field>) -> bool {
    let x_diff_1 = field.sub(&x_vals[1], &x_vals[0]);
    let y_diff_1 = field.sub(&y_vals[1], &y_vals[0]);
    let slope_1 = field.div(&y_diff_1, &x_diff_1);
    let x_diff_2 = field.sub(&x_vals[2], &x_vals[1]);
    let y_diff_2 = field.sub(&y_vals[2], &y_vals[1]);
    let slope_2 = field.div(&y_diff_2, &x_diff_2);
    return slope_1 == slope_2;
  }

  pub fn zeroifier_fft(points: &Vec<u128>, field: &Rc<Field>) -> Polynomial {
    Self::zeroifier_fft_slice(&points[0..], field)
  }

  pub fn zeroifier_fft_slice(points: &[u128], field: &Rc<Field>) -> Polynomial {
    if points.len() == 0 {
      return Polynomial::new(field);
    }
    if points.len() == 1 {
      let mut p = Polynomial::new(field);
      p.term(&field.neg(&points[0]), 0);
      p.term(&1_u128, 1);
      return p;
    }
    let half = points.len() >> 1;
    let left = Self::zeroifier_fft_slice(&points[0..(half)], field);
    let right = Self::zeroifier_fft_slice(&points[(half)..], field);
    return Self::mul_fft(&left, &right, field);
  }

  pub fn zeroifier(points: &Vec<u128>, field: &Rc<Field>) -> Polynomial {
    Self::zeroifier_slice(&points[0..], field)
  }

  pub fn zeroifier_slice(points: &[u128], field: &Rc<Field>) -> Polynomial {
    let mut out = Polynomial::new(field);
    out.term(&1_u128, 0);
    let mut x = Polynomial::new(field);
    x.term(&1_u128, 1);
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
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    poly.term(&f.neg(&9), 0);
    poly.term(&2, 1);
    let mut x_vals: Vec<u128> = Vec::new();
    let mut y_vals: Vec<u128> = Vec::new();
    for i in 0..3 {
      x_vals.push(u128::try_from(i).unwrap());
      y_vals.push(poly.eval(&u128::try_from(i).unwrap()));
    }
    assert!(Polynomial::test_colinearity(&x_vals, &y_vals, &f));
    assert!(Polynomial::test_colinearity_batch(&vec!(x_vals), &vec!(y_vals), &f));
  }

  #[test]
  fn should_compose_polynomial() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let mut root = Polynomial::new(&f);
    root.term(&99, 0);
    root.term(&2, 1);
    root.term(&4, 2);

    let mut inpoly = Polynomial::new(&f);
    inpoly.term(&2, 2);
    inpoly.term(&12, 0);

    let mut expected = Polynomial::new(&f);
    expected.term(&99, 0);
    expected.add(&inpoly.clone().mul_scalar(&2));
    {
      let mut i = inpoly.clone();
      i.exp(2);
      i.mul_scalar(&4);
      expected.add(&i);
    }

    assert!(root.compose(&inpoly).is_equal(&expected));
  }

  #[test]
  fn should_exp_polynomial() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    poly.term(&2, 0);

    for i in 0..10 {
      let mut expected = Polynomial::new(&f);
      expected.term(&f.exp(&2, &u128::try_from(i).unwrap()), 0);
      assert!(poly.clone().exp(i).is_equal(&expected));
    }
  }

  #[test]
  fn should_scale_polynomial() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    poly.term(&1, 0);
    poly.term(&1, 1);
    poly.term(&1, 2);
    poly.term(&1, 3);
    poly.term(&1, 4);

    let mut expected = Polynomial::new(&f);
    expected.term(&1, 0);
    expected.term(&2, 1);
    expected.term(&4, 2);
    expected.term(&8, 3);
    expected.term(&16, 4);

    assert!(expected.is_equal(poly.scale(&2)));
  }

  #[test]
  fn should_interpolate_lagrange() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let size = 128_u128;
    let gen = f.generator(&(size as u32));
    let mut x_vals: Vec<u128> = Vec::new();
    let mut y_vals: Vec<u128> = Vec::new();
    for i in 0..size {
      x_vals.push(f.exp(&gen, &i));
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
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    poly.term(&1, 0);

    let mut divisor = Polynomial::new(&f);
    divisor.term(&0, 3);
    poly.div(&divisor);
  }

  #[test]
  fn should_divide_polynomial() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let mut poly1 = Polynomial::new(&f);
    poly1.term(&1, 2);
    poly1.term(&2, 1);
    poly1.term(&f.neg(&7), 0);

    let mut poly2 = Polynomial::new(&f);
    poly2.term(&1, 1);
    poly2.term(&f.neg(&2), 0);


    let mut expected_q = Polynomial::new(&f);
    expected_q.term(&1, 1);
    expected_q.term(&4, 0);

    let mut expected_r = Polynomial::new(&f);
    expected_r.term(&1, 0);

    let (q, r) = poly1.div(&poly2);

    assert!(q.is_equal(&expected_q));
    assert!(r.is_equal(&expected_r));
  }

  #[test]
  fn should_divide_random_polynomials() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let mut poly1 = Polynomial::new(&f);
    for i in 0..50 {
      poly1.term(&f.random(), i);
    }

    let mut poly2 = Polynomial::new(&f);
    for i in 0..50 {
      poly2.term(&f.random(), i);
    }

    let (q, r) = poly1.div(&poly2);
    let mut out = q.clone();
    out.mul(&poly2);
    out.add(&r);
    assert!(poly1.is_equal(&out));
  }

  #[test]
  fn should_eval_polynomial() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    // 9x^3 - 4x^2 - 20
    let mut poly = Polynomial::new(&f);
    poly.term(&9, 3);
    poly.term(&f.neg(&4), 2);
    poly.term(&f.neg(&20), 0);

    assert_eq!(f.neg(&20), poly.eval(&0));
    assert_eq!(f.neg(&15), poly.eval(&1));
  }

  #[test]
  fn should_eval_polynomial_with_batch_fast() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    for i in 0..50 {
      poly.term(&f.random(), i);
    }

    let size = 2_u128.pow(7);
    let mut g_domain: Vec<u128> = Vec::new();
    for i in 0..size {
      g_domain.push(i);
    }

    let actual = poly.eval_batch(&g_domain);
    let out = poly.eval_batch_fast(&g_domain);
    for i in 0..usize::try_from(size).unwrap() {
      assert_eq!(actual[i], out[i]);
    }
  }

  #[test]
  fn should_eval_polynomial_with_fft() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let mut poly = Polynomial::new(&f);
    for i in 0..50 {
      poly.term(&f.random(), i);
    }

    let size = 2_u128.pow(8);
    let sub_g = f.generator(&(size as u32));
    let mut g_domain: Vec<u128> = Vec::new();
    for i in 0..size {
      g_domain.push(f.exp(&sub_g, &i));
    }

    let actual = poly.eval_batch(&g_domain);
    let out = poly.eval_batch_fft(&g_domain);
    for i in 0..usize::try_from(size).unwrap() {
      assert_eq!(actual[i], out[i]);
    }
  }

  #[test]
  fn should_check_polynomial_equality() {
    let p = 101_u128;
    let g = 0_u128;
    let f = Rc::new(Field::new(p, g));

    let mut poly1 = Polynomial::new(&f);
    poly1.term(&f.neg(&20), 0);
    poly1.term(&2, 2);
    let mut poly2 = Polynomial::new(&f);
    poly2.term(&81, 0);
    poly2.term(&2, 2);
    poly2.term(&0, 5);
    assert!(poly1.is_equal(&poly2));
  }

  #[test]
  fn should_add_polynomials() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    // 2x^2 - 20
    let mut poly1 = Polynomial::new(&f);
    poly1.term(&f.neg(&20), 0);
    poly1.term(&2, 2);

    // 9x^3 - 4x^2 - 20
    let mut poly2 = Polynomial::new(&f);
    poly2.term(&9, 3);
    poly2.term(&f.neg(&4), 2);
    poly2.term(&f.neg(&20), 0);

    // expect 9x^3 - 2x^2 - 40
    let mut out1 = poly1.clone();
    out1.add(&poly2);
    let mut out2 = poly2.clone();
    out2.add(&poly1);

    let mut expected = Polynomial::new(&f);
    expected.term(&9, 3);
    expected.term(&f.neg(&2), 2);
    expected.term(&f.neg(&40), 0);
    assert!(expected.is_equal(&out1));
    assert!(expected.is_equal(&out2));
  }

  #[test]
  fn should_subtract_polynomials() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    // 2x^2 - 20
    let mut poly1 = Polynomial::new(&f);
    poly1.term(&f.neg(&20), 0);
    poly1.term(&2, 2);

    // 9x^3 - 4x^2 - 20
    let mut poly2 = Polynomial::new(&f);
    poly2.term(&9, 3);
    poly2.term(&f.neg(&4), 2);
    poly2.term(&f.neg(&20), 0);

    let mut out1 = poly1.clone();
    out1.sub(&poly2);
    let mut out2 = poly2.clone();
    out2.sub(&poly1);

    let mut expected1 = Polynomial::new(&f);
    expected1.term(&f.neg(&9), 3);
    expected1.term(&6, 2);
    assert!(expected1.is_equal(&out1));

    let mut expected2 = Polynomial::new(&f);
    expected2.term(&9, 3);
    expected2.term(&f.neg(&6), 2);
    assert!(expected2.is_equal(&out2));
  }

  #[test]
  fn should_multiply_polynomials() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    // 2x^2 - 20
    let mut poly1 = Polynomial::new(&f);
    poly1.term(&f.neg(&20), 0);
    poly1.term(&2, 2);

    // 9x^3 - 4x^2 - 20
    let mut poly2 = Polynomial::new(&f);
    poly2.term(&9, 3);
    poly2.term(&f.neg(&4), 2);
    poly2.term(&f.neg(&20), 0);

    let mut poly3 = poly1.clone();
    poly3.mul(&poly2);

    let mut expected = Polynomial::new(&f);
    expected.term(&18, 5);
    expected.term(&f.neg(&8), 4);
    expected.term(&f.neg(&180), 3);
    expected.term(&40, 2);
    expected.term(&400, 0);
    assert!(poly3.is_equal(&expected));
  }

  #[test]
  fn should_multiply_polynomials_fft() {
    let p = 3221225473_u128;
    let g = 5_u128;
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
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let s = 128;
    let size = s as u128;
    let domain_g = f.generator(&(size as u32));
    let mut domain: Vec<u128> = Vec::new();
    for i in 0..size {
      domain.push(f.exp(&domain_g, &i));
    }

    let poly = Polynomial::zeroifier(&domain, &f);
    for i in 0..s {
      assert_eq!(poly.eval(&domain[i]), 0);
    }

    let mut is_zero = true;
    for i in poly.coefs() {
      if i != &0 {
        is_zero = false;
      }
    }
    assert!(!is_zero);
  }

  #[test]
  fn should_build_zeroifier_polynomial_domain() {
    let p = 3221225473_u128;
    let g = 5_u128;
    let f = Rc::new(Field::new(p, g));

    let size = 128_u32;
    let generator = f.generator(&size);
    let domain = f.domain(&generator, size);
    let zeroifier = Polynomial::zeroifier(&domain, &f);
    let mut zeroifier_fft = Polynomial::zeroifier_fft(&domain, &f);
    zeroifier_fft.trim();

    for v in domain {
      assert_eq!(0, zeroifier_fft.eval(&v));
      assert_eq!(0, zeroifier.eval(&v));
    }

    assert!(zeroifier.is_equal(&zeroifier_fft));
  }
}
