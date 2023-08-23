use num_bigint::{BigInt, Sign, BigUint};
use num_bigint::{ToBigInt};
use crate::field::Field;
use std::rc::Rc;
use crate::polynomial::Polynomial;
use crate::channel::Channel;
use crate::mpolynomial::MPolynomial;
use crate::tree::Tree;
use crate::fri::{Fri, FriOptions};
use std::cmp::Ordering;
use std::collections::HashMap;

pub struct Stark {
  offset: BigInt,
  field: Rc<Field>,
  randomizer_count: u32,
  register_count: u32,
  original_trace_len: u32,
  randomized_trace_len: u32,
  expansion_factor: u32,
  omega: BigInt,
  omega_domain: Vec<BigInt>,
  fri_domain_len: u32,
  omicron: BigInt,
  omicron_domain_len: u32,
  omicron_domain: Vec<BigInt>,
  fri: Fri
}

impl Stark {
  pub fn new(
    offset: &BigInt,
    field: &Rc<Field>,
    register_count: u32,
    original_trace_len: u32,
    expansion_factor: u32,
    colinearity_test_count: u32,
    transition_constraints_degree: u32
  ) -> Stark {
    let randomizer_count = 4 * colinearity_test_count;
    let trace_bits = BigInt::from((original_trace_len + randomizer_count)*transition_constraints_degree).bits();
    let omicron_domain_len = 2_u32.pow(u32::try_from(trace_bits).unwrap());
    let fri_domain_len = omicron_domain_len * expansion_factor;
    let omega = field.generator(&BigInt::from(fri_domain_len));
    let omicron = field.generator(&BigInt::from(omicron_domain_len));

    let fri = Fri::new(&FriOptions {
      offset: offset.clone(),
      omega: omega.clone(),
      domain_len: fri_domain_len,
      expansion_factor,
      colinearity_test_count,
    }, field);

    Stark {
      offset: offset.clone(),
      field: Rc::clone(field),
      randomizer_count,
      register_count,
      original_trace_len,
      randomized_trace_len: original_trace_len + randomizer_count,
      expansion_factor,
      omicron_domain_len,
      omicron: omicron.clone(),
      omega_domain: field.domain(&omega, fri_domain_len),
      omega,
      omicron_domain: field.domain(&omicron, omicron_domain_len),
      fri,
      fri_domain_len
    }
  }

  pub fn field(&self) -> &Rc<Field> {
    &self.field
  }

  fn transition_degree_bounds(&self, constraints: &Vec<MPolynomial>) -> Vec<u32> {
    let degree: u32 = self.original_trace_len + self.randomizer_count - 1;
    let mut point_degrees: Vec<u32> = vec!(degree; usize::try_from(1+2*self.register_count).unwrap());
    point_degrees[0] = 1;
    let mut out: Vec<u32> = Vec::new();
    for t in constraints {
      let mut largest = 0;
      for (exps, _) in t.exps() {
        let mut sum = 0;
        for i in 0..point_degrees.len() {
          sum += exps.get(i).unwrap_or(&0) * point_degrees[i];
        }
        if sum > largest {
          largest = sum;
        }
      }
      out.push(largest);
    }
    out
  }

  fn transition_quotient_degree_bounds(&self, constraints: &Vec<MPolynomial>) -> Vec<u32> {
    self.transition_degree_bounds(constraints).iter().map(|v| v - (self.original_trace_len - 1)).collect()
  }

  fn max_degree(&self, constraints: &Vec<MPolynomial>) -> u32 {
    let mut degree: u32 = 0;
    for v in self.transition_quotient_degree_bounds(constraints) {
      degree = std::cmp::max(v, degree);
    }
    let mut max_degree = degree.ilog2();
    if 2_u32.pow(max_degree) != degree {
      max_degree += 1;
    }
    return 2_u32.pow(max_degree) - 1;
  }

  fn transition_zeroifier(&self) -> Polynomial {
    let points = &self.omicron_domain[0..usize::try_from(self.original_trace_len-1).unwrap()];
    Polynomial::zeroifier_fft_slice(points, &self.field)
  }

  fn boundary_zeroifiers(&self, boundary: &Vec<(u32, u32, BigInt)>) -> Vec<Polynomial> {
    let mut zeroifiers: Vec<Polynomial> = Vec::new();
    for i in 0..self.register_count {
      let points: Vec<BigInt> = boundary
        .iter()
        .filter(|(_c, r, _v)| r == &i)
        .map(|(c, _r, _v)| self.omicron_domain[usize::try_from(c.clone()).unwrap()].clone())
        .collect();
      zeroifiers.push(Polynomial::zeroifier_fft(&points, &self.field));
    }
    zeroifiers
  }

  fn boundary_interpolants(&self, boundary: &Vec<(u32, u32, BigInt)>) -> Vec<Polynomial> {
    let mut interpolants: Vec<Polynomial> = Vec::new();
    for i in 0..self.register_count {
      let mut domain: Vec<BigInt> = Vec::new();
      let mut values: Vec<BigInt> = Vec::new();
      for (_, (c, r, v)) in boundary.iter().enumerate() {
        if r != &u32::try_from(i).unwrap() {
          continue;
        }
        domain.push(self.omicron_domain[usize::try_from(c.clone()).unwrap()].clone());
        // domain.push(self.field.exp(&self.omicron, &BigInt::from(c.clone())));
        values.push(v.clone());
      }
      // interpolants.push(Polynomial::lagrange(&domain, &values, &self.field));
      interpolants.push(Polynomial::interpolate_fft(&domain, &values, &self.field));
    }
    interpolants
  }

  fn boundary_quotient_degree_bounds(&self, random_trace_len: u32, boundary: &Vec<(u32, u32, BigInt)>) -> Vec<u32> {
    let random_trace_degree = random_trace_len - 1;
    self.boundary_zeroifiers(boundary).iter().map(|p| random_trace_degree - u32::try_from(p.degree()).unwrap()).collect()
  }

  fn sample_weights(&self, count: u32, randomness: &BigInt) -> Vec<BigInt> {
    let mut out: Vec<BigInt> = Vec::new();
    for i in 0..count {
      let mut hasher = blake3::Hasher::new();
      hasher.update(&randomness.to_bytes_le().1);
      hasher.update(&BigUint::from(i).to_bytes_le());
      let v = BigInt::from_bytes_le(Sign::Plus, hasher.finalize().as_bytes());
      out.push(self.field.modd(&v));
    }
    out
  }

  pub fn prove(&self, trace: &Vec<Vec<BigInt>>, transition_constraints: &Vec<MPolynomial>, boundary: &Vec<(u32, u32, BigInt)>) -> String {
    let mut trace = trace.clone();
    let mut channel = Channel::new();

    for _ in 0..self.randomizer_count {
      let mut r: Vec<BigInt> = Vec::new();
      for _ in 0..self.register_count {
        r.push(self.field.random());
      }
      trace.push(r);
    }

    let mut trace_domain: Vec<BigInt> = vec!(BigInt::from(0); trace.len());
    trace_domain.clone_from_slice(&self.omicron_domain[0..trace.len()]);

    let mut y_vals: Vec<Vec<BigInt>> = Vec::new();
    for i in 0..self.register_count {
      let trace_vals: Vec<BigInt> = trace
        .iter()
        .map(|registers| registers[usize::try_from(i).unwrap()].clone())
        .collect();
      y_vals.push(trace_vals);
    }
    let trace_polys = Polynomial::interpolate_fft_batch(&trace_domain, &y_vals[0..], &self.field);

    let boundary_interpolants = self.boundary_interpolants(&boundary);
    let boundary_zeroifiers = self.boundary_zeroifiers(&boundary);
    let mut boundary_quotients: Vec<Polynomial> = Vec::new();
    for i in 0..usize::try_from(self.register_count).unwrap() {
      let interpolant = &boundary_interpolants[i];
      let zeroifier = &boundary_zeroifiers[i];
      let mut q = trace_polys[i].clone();
      q.sub(interpolant);
      boundary_quotients.push(Polynomial::div_coset(&q, &zeroifier, &self.offset, &self.omega, self.fri_domain_len, &self.field))
      // boundary_quotients.push(q.safe_div(zeroifier));
    }

    let mut boundary_quotient_codewords: Vec<Vec<BigUint>> = Vec::new();
    let mut boundary_quotient_trees: Vec<Tree> = Vec::new();
    for i in 0..usize::try_from(self.register_count).unwrap() {
      let codewords: Vec<BigUint> = boundary_quotients[i].eval_batch_coset(&self.fri.offset, self.fri_domain_len)
        .iter()
        .map(|v| v.to_biguint().unwrap())
        .collect();
      let tree = Tree::build(&codewords);
      channel.push_single(&tree.root());
      boundary_quotient_codewords.push(codewords);
      boundary_quotient_trees.push(tree);
    }

    let mut p_x = Polynomial::new(&self.field);
    p_x.term(&BigInt::from(1), 1);
    let p_x = p_x;

    let mut point: Vec<Polynomial> = Vec::new();
    point.push(p_x.clone());
    point.extend(trace_polys.clone());
    point.extend(trace_polys.iter().map(|p| {
      let mut pp = p.clone();
      pp.scale_precalc(&self.omicron, &self.omicron_domain);
      pp
    }));

    let transition_polynomials: Vec<Polynomial> = transition_constraints.iter().map(|p| p.eval_symbolic(&point)).collect();
    let transition_zeroifier: Polynomial = self.transition_zeroifier();
    let transition_quotients: Vec<Polynomial> = transition_polynomials.iter().map(|p| {
      // p.safe_div(&transition_zeroifier)
      Polynomial::div_coset(&p, &transition_zeroifier, &self.offset, &self.omega, self.fri_domain_len, &self.field)
    }).collect();

    let mut randomizer_poly = Polynomial::new(&self.field);
    let transition_max_degree = self.max_degree(&transition_constraints);
    for i in 0..(transition_max_degree+1) {
      randomizer_poly.term(&self.field.random(), i);
    }

    let randomizer_codeword = randomizer_poly.eval_batch_coset(&self.fri.offset, self.fri_domain_len)
      .iter()
      .map(|v| v.to_biguint().unwrap())
      .collect();
    let randomizer_tree = Tree::build(&randomizer_codeword);
    let randomizer_root = randomizer_tree.root();
    channel.push_single(&randomizer_root);

    let count = u32::try_from(1 + 2 * transition_quotients.len() + 2 * boundary_quotients.len()).unwrap();
    let weights = self.sample_weights(count, &channel.prover_hash().to_bigint().unwrap());

    let bounds = self.transition_quotient_degree_bounds(&transition_constraints);
    for i in 0..bounds.len() {
      if transition_quotients[i].degree() != usize::try_from(bounds[i]).unwrap() {
        panic!("transition quotient degrees do not match expected value");
      }
    }

    let transition_quotient_degree_bounds = self.transition_quotient_degree_bounds(&transition_constraints);
    let boundary_quotient_degree_bounds = self.boundary_quotient_degree_bounds(u32::try_from(trace.len()).unwrap(), &boundary);

    let mut terms: Vec<Polynomial> = Vec::new();
    terms.push(randomizer_poly);
    for i in 0..transition_quotients.len() {
      terms.push(transition_quotients[i].clone());
      let shift = transition_max_degree - transition_quotient_degree_bounds[i];
      let mut shifted = Polynomial::new(&self.field);
      shifted.term(&BigInt::from(1), shift);
      shifted.mul(&transition_quotients[i]);
      terms.push(shifted);
    }
    for i in 0..(usize::try_from(self.register_count).unwrap()) {
      terms.push(boundary_quotients[i].clone());
      let shift = transition_max_degree - boundary_quotient_degree_bounds[i];
      let mut shifted = Polynomial::new(&self.field);
      shifted.term(&BigInt::from(1), shift);
      shifted.mul(&boundary_quotients[i]);
      terms.push(shifted);
    }

    let mut combination = Polynomial::new(&self.field);
    for i in 0..weights.len() {
      let mut w_poly = Polynomial::new(&self.field);
      w_poly.term(&weights[i], 0);
      w_poly.mul(&terms[i]);
      combination.add(&w_poly);
    }

    let combined_codeword = combination.eval_batch_coset(&self.fri.offset, self.fri_domain_len)
      .iter()
      .map(|v| v.to_biguint().unwrap())
      .collect();
    let mut indices = self.fri.prove(&combined_codeword, & mut channel);
    indices.sort_by(|a, b| {
      if a > b {
        return Ordering::Greater;
      }
      return Ordering::Less;
    });
    let mut duplicated_indices = indices.clone();
    duplicated_indices.extend(indices.iter().map(|v| (v + self.expansion_factor) % self.fri_domain_len).collect::<Vec<u32>>());
    let mut quadrupled_indices = duplicated_indices.clone();
    quadrupled_indices.extend(duplicated_indices.iter().map(|v| (v + self.fri_domain_len/2) % self.fri_domain_len).collect::<Vec<u32>>());
    quadrupled_indices.sort_by(|a, b| {
      if a > b {
        return Ordering::Greater;
      }
      return Ordering::Less;
    });

    for (i, bqc) in boundary_quotient_codewords.iter().enumerate() {
      for index in quadrupled_indices.clone() {
        channel.push_single(&bqc[usize::try_from(index).unwrap()]);
        let (path, _) = boundary_quotient_trees[i].open(index);
        channel.push(&path);
      }
    }

    for index in quadrupled_indices {
      channel.push_single(&randomizer_codeword[usize::try_from(index).unwrap()]);
      let (path, _) = randomizer_tree.open(index);
      channel.push(&path);
    }

    channel.serialize()
  }

  pub fn verify(&self, proof: &String, transition_constraints: &Vec<MPolynomial>, boundary: &Vec<(u32, u32, BigInt)>) -> bool {
    let mut channel = Channel::deserialize(proof);
    let mut original_trace_len = 0;
    for (c, _, _) in boundary {
      if c > &original_trace_len {
        original_trace_len = c.clone();
      }
    }
    original_trace_len += 1;

    let randomized_trace_len = original_trace_len + self.randomizer_count;

    let mut boundary_quotient_roots: Vec<BigUint> = Vec::new();
    for _ in 0..self.register_count {
      boundary_quotient_roots.push(channel.pull().data[0].clone());
    }

    let randomizer_root = channel.pull().data[0].clone();

    let count = u32::try_from(1 + 2 * transition_constraints.len() + 2 * usize::try_from(self.register_count).unwrap()).unwrap();
    let weights = self.sample_weights(count, &channel.verifier_hash().to_bigint().unwrap());

    let mut polynomial_vals = self.fri.verify(& mut channel);
    polynomial_vals.sort_by(|(ax, _ay), (bx, _by)| {
      if ax > bx {
        return Ordering::Greater;
      }
      return Ordering::Less;
    });

    let indices: Vec<u32> = polynomial_vals.iter().map(|(x, _y)| x.clone()).collect();
    let values: Vec<BigUint> = polynomial_vals.iter().map(|(_x, y)| y.clone()).collect();

    let mut duplicated_indices = indices.clone();
    duplicated_indices.extend(indices.iter().map(|v| (v + self.expansion_factor) % self.fri_domain_len).collect::<Vec<u32>>());
    duplicated_indices.sort_by(|a, b| {
      if a > b {
        return Ordering::Greater;
      }
      return Ordering::Less;
    });

    let mut leaves: Vec<HashMap<u32, BigUint>> = Vec::new();
    for i in 0..self.register_count {
      let mut leaf_map: HashMap<u32, BigUint> = HashMap::new();
      for j in duplicated_indices.clone() {
        leaf_map.insert(j, channel.pull().data[0].clone());
        let path = &channel.pull().data;
        Tree::verify(
          &boundary_quotient_roots[usize::try_from(i).unwrap()],
          j,
          path,
          leaf_map.get(&j).unwrap()
        );
      }
      leaves.push(leaf_map);
    }

    let mut randomizer_map: HashMap<u32, BigUint> = HashMap::new();
    for i in duplicated_indices {
      let val = channel.pull().data[0].clone();
      let path = &channel.pull().data;
      Tree::verify(&randomizer_root, i, path, &val);
      randomizer_map.insert(i, val);
    }

    let transition_quotient_degree_bounds = self.transition_quotient_degree_bounds(&transition_constraints);
    let boundary_quotient_degree_bounds = self.boundary_quotient_degree_bounds(randomized_trace_len, &boundary);
    let boundary_zeroifiers = self.boundary_zeroifiers(&boundary);
    let boundary_interpolants = self.boundary_interpolants(&boundary);
    let transition_zeroifier = self.transition_zeroifier();
    let transition_constraints_max_degree = self.max_degree(&transition_constraints);

    for i in 0..indices.len() {
      let current_index = indices[i];
      let domain_current_index = self.field.mul(&self.field.g(), &self.omega_domain[usize::try_from(current_index).unwrap()]);
      let next_index = (current_index + self.expansion_factor) % self.fri_domain_len;
      let domain_next_index = self.field.mul(&self.field.g(), &self.omega_domain[usize::try_from(next_index).unwrap()]);
      let mut current_trace = vec!(BigInt::from(0); usize::try_from(self.register_count).unwrap());
      let mut next_trace = vec!(BigInt::from(0); usize::try_from(self.register_count).unwrap());

      for j in 0..usize::try_from(self.register_count).unwrap() {
        let zeroifier = &boundary_zeroifiers[j];
        let interpolant = &boundary_interpolants[j];

        current_trace[j] = self.field.add(
          &self.field.mul(&leaves[j].get(&current_index).unwrap().to_bigint().unwrap(), &zeroifier.eval(&domain_current_index)),
          &interpolant.eval(&domain_current_index)
        );
        next_trace[j] = self.field.add(
          &self.field.mul(&leaves[j].get(&next_index).unwrap().to_bigint().unwrap(), &zeroifier.eval(&domain_next_index)),
          &interpolant.eval(&domain_next_index)
        );
      }

      let mut point: Vec<BigInt> = Vec::new();
      point.push(domain_current_index.clone());
      point.extend(current_trace.clone());
      point.extend(next_trace.clone());

      let transition_constraint_values: Vec<BigInt> = transition_constraints.iter().map(|c| c.eval(&point)).collect();

      let mut terms: Vec<BigInt> = Vec::new();
      terms.push(randomizer_map.get(&current_index).unwrap().to_bigint().unwrap());
      for j in 0..transition_constraint_values.len() {
        let tcv = &transition_constraint_values[j];
        let q = self.field.div(tcv, &transition_zeroifier.eval(&domain_current_index));
        terms.push(q.clone());
        let shift = transition_constraints_max_degree - transition_quotient_degree_bounds[j];
        terms.push(self.field.mul(&q, &self.field.exp(&domain_current_index, &BigInt::from(shift))));
      }

      for j in 0..usize::try_from(self.register_count).unwrap() {
        let bqv = leaves[j].get(&current_index).unwrap().to_bigint().unwrap();
        terms.push(bqv.clone());
        let shift = transition_constraints_max_degree - boundary_quotient_degree_bounds[j];
        terms.push(self.field.mul(&bqv, &self.field.exp(&domain_current_index, &BigInt::from(shift))));
      }

      let mut combination = BigInt::from(0);
      for j in 0..terms.len() {
        combination = self.field.add(&combination, &self.field.mul(&terms[j], &weights[j]));
      }
      if combination.to_biguint().unwrap() != values[i] {
        panic!("invalid combination value");
      }
    }

    return true;
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn should_make_verify_stark_proof() {
    let p = BigInt::from(1) + BigInt::from(407) * BigInt::from(2).pow(119);
    let g = BigInt::from(85408008396924667383611388730472331217_u128);
    let f = Rc::new(Field::new(p, g.clone()));

    let sequence_len = 40;
    let stark = Stark::new(
      &g.clone(),
      &f,
      2,
      sequence_len,
      128,
      18,
      2
    );

    let mut trace: Vec<Vec<BigInt>> = Vec::new();
    trace.push(vec!(BigInt::from(2), BigInt::from(3)));
    trace.push(vec!(BigInt::from(4), BigInt::from(9)));
    while trace.len() < sequence_len.try_into().unwrap() {
      let e1 = &trace[trace.len() - 1][0];
      let e2 = &trace[trace.len() - 1][1];
      trace.push(vec!(f.mul(e1, e1), f.mul(e2, e2)));
    }

    let boundary_constraints = vec!(
      (0, 0, BigInt::from(2)),
      (0, 1, BigInt::from(3)),
      (sequence_len-1, 0, trace[trace.len()-1][0].clone()),
      (sequence_len-1, 1, trace[trace.len()-1][1].clone())
    );

    let variables = MPolynomial::variables(1+2*2, &f);

    let cycle_index = &variables[0];
    let prev_state = &variables[1..3];
    let next_state = &variables[3..];
    let mut transition_constraints: Vec<MPolynomial> = Vec::new();
    {
      let mut c = prev_state[0].clone();
      c.mul(&prev_state[0]);
      c.sub(&next_state[0]);
      transition_constraints.push(c);
    }
    {
      let mut c = prev_state[1].clone();
      c.mul(&prev_state[1]);
      c.sub(&next_state[1]);
      transition_constraints.push(c);
    }

    let proof = stark.prove(&trace, &transition_constraints, &boundary_constraints);
    stark.verify(&proof, &transition_constraints, &boundary_constraints);
  }
}
