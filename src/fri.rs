use num_bigint::{Sign, BigUint, BigInt};
use num_bigint::ToBigInt;
use crate::field::Field;
use std::rc::Rc;
use crate::tree::Tree;
use crate::channel::Channel;
use crate::polynomial::Polynomial;
use std::collections::HashMap;

pub struct FriOptions {
  pub offset: u128,
  pub omega: u128,
  pub domain_len: u32,
  pub expansion_factor: u32,
  pub colinearity_test_count: u32
}

pub struct Fri {
  pub offset: u128,
  pub omega: u128,
  pub domain_len: u32,
  pub field: Rc<Field>,
  pub expansion_factor: u32,
  pub colinearity_test_count: u32,
  domain: Vec<u128>,
  round_count: u32
}

impl Fri {
  pub fn new(options: &FriOptions, field: &Rc<Field>) -> Fri {
    // calculate number of rounds
    let mut codeword_len = options.domain_len;
    let mut round_count = 0;
    while codeword_len > options.expansion_factor && 4 * options.colinearity_test_count < codeword_len {
      codeword_len /= 2;
      round_count+= 1;
    }
    Fri {
      offset: options.offset.clone(),
      omega: options.omega.clone(),
      domain_len: options.domain_len,
      field: Rc::clone(field),
      expansion_factor: options.expansion_factor,
      colinearity_test_count: options.colinearity_test_count,
      domain: field.coset(options.domain_len, &options.offset),
      round_count
    }
  }

  pub fn domain(&self) -> &Vec<u128> {
    &self.domain
  }

  fn round_count(&self) -> u32 {
    self.round_count
  }

  pub fn prove(&self, codeword: &Vec<u128>, channel: & mut Channel) -> Vec<u32> {
    if self.domain_len != u32::try_from(codeword.len()).unwrap() {
      panic!("initial codeword does not match domain len");
    }
    let codewords = self.commit(codeword, channel);
    let top_indices = Self::sample_indices(
      &self.field.sample_bytes(&channel.prover_hash()),
      codewords[1].len().try_into().unwrap(),
      codewords[codewords.len() - 1].len().try_into().unwrap(),
      self.colinearity_test_count
    );
    let mut indices: Vec<u32> = top_indices.clone();
    let codeword_trees: Vec<Tree> = codewords.iter().map(|word| Tree::build(&Tree::u128s_to_bytes(word))).collect();
    for i in 0..(codewords.len() - 1) {
      indices = indices.iter().map(|index| index % u32::try_from(codewords[i].len() >> 1).unwrap()).collect();
      self.query(&codewords[i], &codewords[i+1], &indices, channel, &codeword_trees[i], &codeword_trees[i+1]);
    }
    top_indices
  }

  fn query(
    &self,
    current_codeword: &Vec<u128>,
    next_codeword: &Vec<u128>,
    indices_c: &Vec<u32>,
    channel: & mut Channel,
    current_codeword_tree: &Tree,
    next_codeword_tree: &Tree
  ) {
    let indices_a: Vec<u32> = indices_c.to_vec();
    let indices_b: Vec<u32> = indices_c.iter().map(|val| val + ((current_codeword.len() >> 1) as u32)).collect();
    for i in 0..usize::try_from(self.colinearity_test_count).unwrap() {
      channel.push(&vec!(
        current_codeword[usize::try_from(indices_a[i]).unwrap()].clone(),
        current_codeword[usize::try_from(indices_b[i]).unwrap()].clone(),
        next_codeword[usize::try_from(indices_c[i]).unwrap()].clone()
      ));
    }
    for i in 0..usize::try_from(self.colinearity_test_count).unwrap() {
      channel.push_path(&current_codeword_tree.open(indices_a[i]).0);
      channel.push_path(&current_codeword_tree.open(indices_b[i]).0);
      channel.push_path(&next_codeword_tree.open(indices_c[i]).0);
    }
  }

  fn commit(&self, codeword: &Vec<u128>, channel: & mut Channel) -> Vec<Vec<u128>> {
    let mut codewords: Vec<Vec<u128>> = Vec::new();
    let mut codeword = codeword.clone();
    let two_inv = self.field.inv(&2);

    // invert the entire domain using repeated multiplications
    // e.g. 1/4 = (1/2) * (1/2)
    // 1/x^2 = (1/x) * (1/x)
    let inv_offset = self.field.inv(&self.offset);
    let inv_offset_domain = self.field.domain(&inv_offset, 2_u32.pow(self.round_count()));

    let inv_omega = self.field.inv(&self.omega);
    let inv_domain = self.field.domain(&inv_omega, self.domain_len);

    let mut exp: usize = 1;

    for x in 0..self.round_count() {
      let root = Tree::commit_elements(&codeword);
      channel.push_root(&root);
      if x == self.round_count() - 1 {
        break;
      }
      codewords.push(codeword.clone());
      // now split the last codeword and fold into a set
      // of points from a polynomial of half the degree
      // of the previous codewords, similar to a FFT
      let alpha = self.field.sample_bytes(&channel.prover_hash());
      let next_len = codeword.len() >> 1;
      codeword = codeword[0..next_len].iter().enumerate().map(|(index, val)| {
        let inv_omega = self.field.mul(&inv_offset_domain[exp], &inv_domain[exp * index]);
        // ( (one + alpha / (offset * (omega^i)) ) * codeword[i]
        let a = self.field.mul(&val, &self.field.add(&Field::one(), &self.field.mul(&alpha, &inv_omega)));
        //  (one - alpha / (offset * (omega^i)) ) * codeword[len(codeword)//2 + i] ) for i in range(len(codeword)//2)]
        let b = self.field.mul(
          &self.field.sub(&Field::one(), &self.field.mul(&alpha, &inv_omega)),
          &codeword[(codeword.len() >> 1) + index]
        );
        return self.field.mul(&two_inv, &self.field.add(&a, &b));
      }).collect();

      exp *= 2;
    }
    channel.push(&codeword);
    codewords.push(codeword);
    codewords
  }

  fn sample_indices(seed: &u128, size: u32, reduced_size: u32, count: u32) -> Vec<u32> {
    if count > 2*reduced_size {
      panic!("not enough entropy");
    }
    if count > reduced_size {
      panic!("cannot sample more indices than available");
    }

    let mut indices: Vec<u32> = Vec::new();
    let mut reduced_indices: HashMap<u32, bool> = HashMap::new();
    let mut counter: u32 = 0;
    let f = Field::new(u128::try_from(size).unwrap(), 0);
    while indices.len() < (count as usize) {
      let mut hasher = blake3::Hasher::new();
      hasher.update(&seed.to_le_bytes());
      hasher.update(&counter.to_le_bytes());
      let v = f.sample_bytes(hasher.finalize().as_bytes());
      // let v = BigInt::from_bytes_le(Sign::Plus, hasher.finalize().as_bytes());
      let index = Field::bigint_to_u32(&v);
      let reduced_index = index % reduced_size;
      counter += 1;
      if !reduced_indices.contains_key(&reduced_index) {
        indices.push(index);
        reduced_indices.insert(reduced_index, true);
      }
    }
    indices
  }

  pub fn verify(&self, channel: & mut Channel) -> Vec<(u32, u128)> {
    let mut out: Vec<(u32, u128)> = Vec::new();
    let mut omega = self.omega.clone();
    let mut offset = self.offset.clone();

    let mut roots: Vec<[u8; 32]> = Vec::new();
    let mut alphas: Vec<u128> = Vec::new();

    for _ in 0..self.round_count() {
      roots.push(channel.pull_root().clone());
      alphas.push(self.field.sample_bytes(&channel.verifier_hash()));
    }

    let last_codeword = channel.pull_u128s().clone();
    if roots[roots.len() - 1] != Tree::commit_elements(&last_codeword) {
      panic!("last codeword root mismatch");
    }

    let degree: usize = (last_codeword.len() / usize::try_from(self.expansion_factor).unwrap()) - 1;
    let mut last_omega = omega.clone();
    let mut last_offset = offset.clone();
    for _ in 0..(self.round_count() - 1) {
      last_omega = self.field.mul(&last_omega, &last_omega);
      last_offset = self.field.mul(&last_offset, &last_offset);
    }
    if self.field.inv(&last_omega) != self.field.exp(&last_omega, &u128::try_from(last_codeword.len() - 1).unwrap()) {
      panic!("omega order incorrect");
    }

    let last_domain: Vec<u128> = last_codeword.iter().enumerate().map(|(index, _)| {
      return self.field.mul(&last_offset, &self.field.exp(&last_omega, &u128::try_from(index).unwrap()));
    }).collect();

    let poly = Polynomial::lagrange(&last_domain, &last_codeword, &self.field);
    for i in 0..last_domain.len() {
      if poly.eval(&last_domain[i]) != last_codeword[i] {
        println!("{} {}", poly.eval(&last_domain[i]), last_codeword[i]);
        panic!("interpolated polynomial is incorrect");
      }
    }
    if poly.degree() > degree {
      panic!("last codeword does not match polynomial of low enough degree");
    }

    let top_indices = Self::sample_indices(
      &self.field.sample_bytes(&channel.verifier_hash()),
      self.domain_len >> 1,
      self.domain_len >> (self.round_count() - 1),
      self.colinearity_test_count
    );
    let mut colinearity_x_vals: Vec<Vec<u128>> = Vec::new();
    let mut colinearity_y_vals: Vec<Vec<u128>> = Vec::new();
    for i in 0..usize::try_from(self.round_count() - 1).unwrap() {
      let indices_c: Vec<u32> = top_indices.iter().map(|val| val % (self.domain_len >> (i+1))).collect();
      let indices_a: Vec<u32> = indices_c.clone();
      let indices_b: Vec<u32> = indices_a.iter().map(|val| val + (self.domain_len >> (i+1))).collect();

      let mut aa: Vec<[u8; 32]> = Vec::new();
      let mut bb: Vec<[u8; 32]> = Vec::new();
      let mut cc: Vec<[u8; 32]> = Vec::new();
      for j in 0..usize::try_from(self.colinearity_test_count).unwrap() {
        let y_points_msg = channel.pull_u128s();
        let ay = y_points_msg[0].clone();
        let by = y_points_msg[1].clone();
        let cy = y_points_msg[2].clone();
        aa.push(Tree::u128_to_bytes(&ay));
        bb.push(Tree::u128_to_bytes(&by));
        cc.push(Tree::u128_to_bytes(&cy));
        if i == 0 {
          out.push((indices_a[j], y_points_msg[0].clone()));
          out.push((indices_b[j], y_points_msg[1].clone()));
        }

        let ax = self.field.mul(&offset, &self.field.exp(&omega, &u128::try_from(indices_a[j]).unwrap()));
        let bx = self.field.mul(&offset, &self.field.exp(&omega, &u128::try_from(indices_b[j]).unwrap()));

        let cx = alphas[usize::try_from(i).unwrap()].clone();

        colinearity_x_vals.push(vec!(ax, bx, cx));
        colinearity_y_vals.push(vec!(ay, by, cy));
      }

      if !Polynomial::test_colinearity_batch(&colinearity_x_vals, &colinearity_y_vals, &self.field) {
        panic!("colinearity test failed");
      }

      for j in 0..usize::try_from(self.colinearity_test_count).unwrap() {
        Tree::verify(&roots[i], indices_a[j], &channel.pull_path(), &aa[j]);
        Tree::verify(&roots[i], indices_b[j], &channel.pull_path(), &bb[j]);
        Tree::verify(&roots[i+1], indices_c[j], &channel.pull_path(), &cc[j]);
      }

      omega = self.field.mul(&omega, &omega);
      offset = self.field.mul(&offset, &offset);
    }

    out
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn should_make_verify_fri_proof() {
    let mut channel = Channel::new();
    let p = 1 + 407 * 2_u128.pow(119);
    let g = 85408008396924667383611388730472331217_u128;
    // let p = 3221225473_u128;
    // let g = 5_u128;
    let f = Rc::new(Field::new(p, g.clone()));
    let domain_size: u32 = 64;
    let domain_g = f.generator(&domain_size);

    let fri = Fri::new(&FriOptions {
      offset: g.clone(),
      omega: domain_g.clone(),
      domain_len: domain_size,
      expansion_factor: 4,
      colinearity_test_count: 4
    }, &f);

    let mut poly = Polynomial::new(&f);
    poly.term(&3, 2);
    let mut points: Vec<u128> = Vec::new();
    for i in fri.domain() {
      points.push(poly.eval(&i));
    }
    fri.prove(&points, & mut channel);
    fri.verify(& mut channel);
  }
}
