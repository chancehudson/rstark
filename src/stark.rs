use crate::channel::Channel;
use crate::field::Field;
use crate::fri::{Fri, FriOptions};
use crate::mpolynomial::MPolynomial;
use crate::polynomial::Polynomial;
use crate::tree::Tree;
use crate::FieldElement;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::rc::Rc;

pub struct Stark<T: FieldElement> {
    offset: T,
    field: Rc<Field<T>>,
    randomizer_count: u32,
    register_count: u32,
    original_trace_len: u32,
    expansion_factor: u32,
    omega: T,
    omega_domain: Vec<T>,
    fri_domain_len: u32,
    omicron: T,
    omicron_domain: Vec<T>,
    fri: Fri<T>,
}

impl<T: FieldElement> Stark<T> {
    pub fn new(
        offset: &T,
        field: &Rc<Field<T>>,
        register_count: u32,
        original_trace_len: u32,
        expansion_factor: u32,
        colinearity_test_count: u32,
        transition_constraints_degree: u32,
    ) -> Stark<T> {
        let randomizer_count = 4 * colinearity_test_count;
        let trace_bits = T::from_u32(
            (original_trace_len + randomizer_count) * transition_constraints_degree,
            field.p(),
        )
        .bits();
        let omicron_domain_len = 2_u32.pow(u32::try_from(trace_bits).unwrap());
        let fri_domain_len = omicron_domain_len * expansion_factor;
        let (omega, _) = field.generator_cache(&fri_domain_len);
        let (omicron, _) = field.generator_cache(&omicron_domain_len);

        let fri = Fri::new(
            &FriOptions {
                offset: offset.clone(),
                omega: omega.clone(),
                domain_len: fri_domain_len,
                expansion_factor,
                colinearity_test_count,
            },
            field,
        );

        Stark {
            offset: offset.clone(),
            field: Rc::clone(field),
            randomizer_count,
            register_count,
            original_trace_len,
            expansion_factor,
            omicron: omicron.clone(),
            omega_domain: field.domain(&omega, fri_domain_len),
            omega,
            omicron_domain: field.domain(&omicron, omicron_domain_len),
            fri,
            fri_domain_len,
        }
    }

    pub fn field(&self) -> &Rc<Field<T>> {
        &self.field
    }

    fn transition_degree_bounds(&self, constraint: &MPolynomial<T>) -> u32 {
        let degree: u32 = self.original_trace_len + self.randomizer_count - 1;
        let mut point_degrees: Vec<u32> =
            vec![degree; usize::try_from(1 + 2 * self.register_count).unwrap()];
        point_degrees[0] = 1;
        let mut out: u32 = 0;
        for (exps, _) in constraint.exps() {
            let mut sum = 0;
            for (i, v) in point_degrees.iter().enumerate() {
                sum += exps.get(i).unwrap_or(&0) * v;
            }
            if sum > out {
                out = sum;
            }
        }
        out
    }

    fn transition_quotient_degree_bound(&self, constraint: &MPolynomial<T>) -> u32 {
        self.transition_degree_bounds(constraint) - (self.original_trace_len - 1)
    }

    fn max_degree(&self, constraint: &MPolynomial<T>) -> u32 {
        let degree: u32 = self.transition_quotient_degree_bound(constraint);
        let mut max_degree = degree.ilog2();
        if 2_u32.pow(max_degree) != degree {
            max_degree += 1;
        }
        2_u32.pow(max_degree) - 1
    }

    fn transition_zeroifier(&self) -> Polynomial<T> {
        let points = &self.omicron_domain[0..usize::try_from(self.original_trace_len - 1).unwrap()];
        Polynomial::zeroifier_fft_slice(points, &self.field)
    }

    fn boundary_zeroifiers(&self, boundary: &Vec<(u32, u32, T)>) -> Vec<Polynomial<T>> {
        let mut zeroifiers = Vec::new();
        for i in 0..self.register_count {
            let points: Vec<T> = boundary
                .iter()
                .filter(|(_c, r, _v)| r == &i)
                .map(|(c, _r, _v)| self.omicron_domain[usize::try_from(*c).unwrap()].clone())
                .collect();
            if points.is_empty() {
                let mut p = Polynomial::new(&self.field);
                p.term(&self.omicron, 0);
                zeroifiers.push(p)
            } else {
                zeroifiers.push(Polynomial::zeroifier_fft(&points, &self.field));
            }
        }
        zeroifiers
    }

    fn boundary_interpolants(&self, boundary: &Vec<(u32, u32, T)>) -> Vec<Polynomial<T>> {
        let mut interpolants: Vec<Polynomial<T>> = Vec::new();
        for i in 0..self.register_count {
            let mut domain = Vec::new();
            let mut values = Vec::new();
            for (_, (c, r, v)) in boundary.iter().enumerate() {
                if *r != i {
                    continue;
                }
                domain.push(self.omicron_domain[usize::try_from(*c).unwrap()].clone());
                // domain.push(self.field.exp(&self.omicron, &BigInt::from(c.clone())));
                values.push(v.clone());
            }
            // interpolants.push(Polynomial::lagrange(&domain, &values, &self.field));
            interpolants.push(Polynomial::interpolate_fft(&domain, &values, &self.field));
        }
        interpolants
    }

    fn boundary_quotient_degree_bounds(
        &self,
        random_trace_len: u32,
        boundary: &Vec<(u32, u32, T)>,
    ) -> Vec<u32> {
        let random_trace_degree = random_trace_len - 1;
        self.boundary_zeroifiers(boundary)
            .iter()
            .map(|p| random_trace_degree - u32::try_from(p.degree()).unwrap())
            .collect()
    }

    fn sample_weights(&self, count: u32, randomness: &T) -> Vec<T> {
        let mut out = Vec::new();
        for i in 0..count {
            let mut hasher = blake3::Hasher::new();
            hasher.update(&randomness.to_bytes_le());
            hasher.update(&T::from_u32(i, self.field().p()).to_bytes_le());
            let v = T::from_bytes_le(hasher.finalize().as_bytes(), self.field().p());
            out.push(v);
        }
        out
    }

    pub fn prove(
        &self,
        trace: &Vec<Vec<T>>,
        transition_constraints: &Vec<MPolynomial<T>>,
        boundary: &Vec<(u32, u32, T)>,
    ) -> String {
        let mut trace = trace.clone();
        let mut channel = Channel::new();

        for _ in 0..self.randomizer_count {
            let mut r: Vec<T> = Vec::new();
            for _ in 0..self.register_count {
                r.push(self.field.random());
            }
            trace.push(r);
        }

        let mut trace_domain = vec![self.field.zero(); trace.len()];
        trace_domain.clone_from_slice(&self.omicron_domain[0..trace.len()]);

        let mut y_vals = Vec::new();
        for i in 0..self.register_count {
            let trace_vals = trace
                .iter()
                .map(|registers| registers[usize::try_from(i).unwrap()].clone())
                .collect();
            y_vals.push(trace_vals);
        }
        let trace_polys =
            Polynomial::interpolate_fft_batch(&trace_domain, &y_vals[0..], &self.field);

        let boundary_interpolants = self.boundary_interpolants(boundary);
        let boundary_zeroifiers = self.boundary_zeroifiers(boundary);
        let mut boundary_quotients = Vec::new();
        for i in 0..usize::try_from(self.register_count).unwrap() {
            let interpolant = &boundary_interpolants[i];
            let zeroifier = &boundary_zeroifiers[i];
            let mut q = trace_polys[i].clone();
            q.sub(interpolant);
            boundary_quotients.push(Polynomial::div_coset(
                &q,
                zeroifier,
                &self.offset,
                &self.omega,
                self.fri_domain_len,
                &self.field,
            ))
            // boundary_quotients.push(q.safe_div(zeroifier));
        }

        let mut boundary_quotient_codewords = Vec::new();
        let mut boundary_quotient_trees: Vec<Tree<T>> = Vec::new();
        let codewords = Polynomial::eval_batch_batch_coset(
            boundary_quotients.clone(),
            &self.fri.offset,
            self.fri_domain_len,
            &self.field,
        );
        for i in 0..usize::try_from(self.register_count).unwrap() {
            let c = codewords[i].iter().map(|v| v.to_bytes_le_sized()).collect();
            let tree = Tree::build(&c);
            channel.push_single(&tree.root());
            boundary_quotient_codewords.push(c);
            boundary_quotient_trees.push(tree);
        }

        let mut p_x = Polynomial::new(&self.field);
        p_x.term(&self.field.one(), 1);
        let p_x = p_x;

        let mut point = Vec::new();
        point.push(p_x.clone());
        point.extend(trace_polys.clone());
        point.extend(trace_polys.iter().map(|p| {
            let mut pp = p.clone();
            pp.scale_precalc(&self.omicron, &self.omicron_domain);
            pp
        }));

        let mut transition_weights = Vec::new();
        {
            let coef = self
                .field
                .sample(T::from_bytes_le(&channel.prover_hash(), self.field().p()));
            transition_weights.push(coef.clone());
            for i in 1..transition_constraints.len() {
                transition_weights.push(self.field.mul(&transition_weights[i - 1], &coef));
            }
        }

        // combine all transition constraints using a random linear combination
        let mut single_transition_constraint = MPolynomial::new(&self.field);
        for (i, t) in transition_constraints.iter().enumerate() {
            single_transition_constraint.add(t.clone().mul_scalar(&transition_weights[i]));
        }

        let transition_polynomial = single_transition_constraint.eval_symbolic(&point);
        let transition_zeroifier = self.transition_zeroifier();
        let transition_quotient = Polynomial::div_coset(
            &transition_polynomial,
            &transition_zeroifier,
            &self.offset,
            &self.omega,
            self.fri_domain_len,
            &self.field,
        );

        let mut randomizer_poly = Polynomial::new(&self.field);
        let transition_max_degree = self.max_degree(&single_transition_constraint);
        for i in 0..(transition_max_degree + 1) {
            randomizer_poly.term(&self.field.random(), i);
        }

        let randomizer_codeword = randomizer_poly
            .eval_batch_coset(&self.fri.offset, self.fri_domain_len)
            .iter()
            .map(|v| v.to_bytes_le_sized())
            .collect::<Vec<[u8; 32]>>();
        let randomizer_tree: Tree<T> = Tree::build(&randomizer_codeword);
        let randomizer_root = randomizer_tree.root();
        channel.push_single(&randomizer_root);

        let count =
            u32::try_from(1 + 2/*transition_quotients.len()*/ + 2 * boundary_quotients.len())
                .unwrap();
        let weights = self.sample_weights(
            count,
            &T::from_bytes_le(&channel.prover_hash(), self.field().p()),
        );

        let bounds = self.transition_quotient_degree_bound(&single_transition_constraint);
        if transition_quotient.degree() != usize::try_from(bounds).unwrap() {
            panic!("transition quotient degrees do not match expected value");
        }

        let transition_quotient_degree_bound =
            self.transition_quotient_degree_bound(&single_transition_constraint);
        let boundary_quotient_degree_bounds =
            self.boundary_quotient_degree_bounds(u32::try_from(trace.len()).unwrap(), boundary);

        let mut terms = Vec::new();
        terms.push(randomizer_poly);

        terms.push(transition_quotient.clone());
        let shift = transition_max_degree - transition_quotient_degree_bound;
        terms.push(transition_quotient.shift_and_clone(shift));

        for i in 0..(usize::try_from(self.register_count).unwrap()) {
            terms.push(boundary_quotients[i].clone());
            let shift = transition_max_degree - boundary_quotient_degree_bounds[i];
            terms.push(boundary_quotients[i].shift_and_clone(shift));
        }

        let mut combination = Polynomial::new(&self.field);
        for i in 0..weights.len() {
            let mut w_poly = Polynomial::new(&self.field);
            w_poly.term(&weights[i], 0);
            w_poly.mul(&terms[i]);
            combination.add(&w_poly);
        }

        let combined_codeword = combination
            .eval_batch_coset(&self.fri.offset, self.fri_domain_len)
            .to_vec();
        let mut indices = self.fri.prove(&combined_codeword, &mut channel);
        indices.sort_by(|a, b| {
            if a > b {
                return Ordering::Greater;
            }
            Ordering::Less
        });
        let mut duplicated_indices = indices.clone();
        duplicated_indices.extend(
            indices
                .iter()
                .map(|v| (v + self.expansion_factor) % self.fri_domain_len)
                .collect::<Vec<u32>>(),
        );
        let mut quadrupled_indices = duplicated_indices.clone();
        quadrupled_indices.extend(
            duplicated_indices
                .iter()
                .map(|v| (v + self.fri_domain_len / 2) % self.fri_domain_len)
                .collect::<Vec<u32>>(),
        );
        quadrupled_indices.sort_by(|a, b| {
            if a > b {
                return Ordering::Greater;
            }
            Ordering::Less
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

    pub fn verify(
        &self,
        proof: &str,
        transition_constraints: &Vec<MPolynomial<T>>,
        boundary: &Vec<(u32, u32, T)>,
    ) -> bool {
        let mut channel = Channel::deserialize(proof);
        let mut original_trace_len = 0;
        for (c, _, _) in boundary {
            if c > &original_trace_len {
                original_trace_len = *c;
            }
        }
        original_trace_len += 1;

        let randomized_trace_len = original_trace_len + self.randomizer_count;

        let mut boundary_quotient_roots = Vec::new();
        for _ in 0..self.register_count {
            boundary_quotient_roots.push(channel.pull_root());
        }

        let mut transition_weights = Vec::new();
        {
            let coef = self
                .field
                .sample(T::from_bytes_le(&channel.verifier_hash(), self.field().p()));
            transition_weights.push(coef.clone());
            for i in 1..transition_constraints.len() {
                transition_weights.push(self.field.mul(&transition_weights[i - 1], &coef));
            }
        }

        // combine all transition constraints using a random linear combination
        let mut single_transition_constraint = MPolynomial::new(&self.field);
        for (i, t) in transition_constraints.iter().enumerate() {
            single_transition_constraint.add(t.clone().mul_scalar(&transition_weights[i]));
        }

        let randomizer_root = channel.pull_root();

        let count = u32::try_from(
            1 + 2 * transition_constraints.len()
                + 2 * usize::try_from(self.register_count).unwrap(),
        )
        .unwrap();
        let weights = self.sample_weights(
            count,
            &T::from_bytes_le(&channel.verifier_hash(), self.field().p()),
        );

        let mut polynomial_vals = self.fri.verify(&mut channel);
        polynomial_vals.sort_by(|(ax, _ay), (bx, _by)| {
            if ax > bx {
                return Ordering::Greater;
            }
            Ordering::Less
        });

        let indices: Vec<u32> = polynomial_vals.iter().map(|(x, _y)| *x).collect();
        let values: Vec<T> = polynomial_vals.iter().map(|(_x, y)| y.clone()).collect();

        let mut duplicated_indices = indices.clone();
        duplicated_indices.extend(
            indices
                .iter()
                .map(|v| (v + self.expansion_factor) % self.fri_domain_len)
                .collect::<Vec<u32>>(),
        );
        duplicated_indices.sort_by(|a, b| {
            if a > b {
                return Ordering::Greater;
            }
            Ordering::Less
        });

        let mut leaves = Vec::new();
        for i in 0..self.register_count {
            let mut leaf_map = HashMap::new();
            for j in duplicated_indices.clone() {
                leaf_map.insert(j, channel.pull_path()[0]);
                let path = &channel.pull_path();
                Tree::<T>::verify(
                    &boundary_quotient_roots[usize::try_from(i).unwrap()],
                    j,
                    path,
                    leaf_map.get(&j).unwrap(),
                );
            }
            leaves.push(leaf_map);
        }

        let mut randomizer_map = HashMap::new();
        for i in duplicated_indices {
            let val = channel.pull().data.clone().try_into().unwrap();
            let path = &channel.pull_path();
            Tree::<T>::verify(&randomizer_root, i, path, &val);
            randomizer_map.insert(i, val);
        }

        let transition_quotient_degree_bound =
            self.transition_quotient_degree_bound(&single_transition_constraint);
        let boundary_quotient_degree_bounds =
            self.boundary_quotient_degree_bounds(randomized_trace_len, boundary);
        let boundary_zeroifiers = self.boundary_zeroifiers(boundary);
        let boundary_interpolants = self.boundary_interpolants(boundary);
        let transition_zeroifier = self.transition_zeroifier();
        let transition_constraints_max_degree = self.max_degree(&single_transition_constraint);

        for i in 0..indices.len() {
            let current_index = indices[i];
            let domain_current_index = self.field.mul(
                self.field.g(),
                &self.omega_domain[usize::try_from(current_index).unwrap()],
            );
            let next_index = (current_index + self.expansion_factor) % self.fri_domain_len;
            let domain_next_index = self.field.mul(
                self.field.g(),
                &self.omega_domain[usize::try_from(next_index).unwrap()],
            );
            let mut current_trace =
                vec![self.field.zero(); usize::try_from(self.register_count).unwrap()];
            let mut next_trace =
                vec![self.field.zero(); usize::try_from(self.register_count).unwrap()];

            for j in 0..usize::try_from(self.register_count).unwrap() {
                let zeroifier = &boundary_zeroifiers[j];
                let interpolant = &boundary_interpolants[j];

                current_trace[j] = self.field.add(
                    &self.field.mul(
                        &T::from_bytes_le(leaves[j].get(&current_index).unwrap(), self.field().p()),
                        &zeroifier.eval(&domain_current_index),
                    ),
                    &interpolant.eval(&domain_current_index),
                );
                next_trace[j] = self.field.add(
                    &self.field.mul(
                        &T::from_bytes_le(leaves[j].get(&next_index).unwrap(), self.field().p()),
                        &zeroifier.eval(&domain_next_index),
                    ),
                    &interpolant.eval(&domain_next_index),
                );
            }

            let mut point = Vec::new();
            point.push(domain_current_index.clone());
            point.extend(current_trace.clone());
            point.extend(next_trace.clone());

            let transition_constraint_value = single_transition_constraint.eval(&point);
            let transition_zeroifier_eval_inv = self
                .field
                .inv(&transition_zeroifier.eval(&domain_current_index));

            let mut terms = Vec::new();
            terms.push(T::from_bytes_le(
                randomizer_map.get(&current_index).unwrap(),
                self.field().p(),
            ));

            // power map for domain_current_index
            let mut power_map = HashMap::new();

            let q = self
                .field
                .mul(&transition_constraint_value, &transition_zeroifier_eval_inv);
            terms.push(q.clone());
            let shift = transition_constraints_max_degree - transition_quotient_degree_bound;
            {
                let exp = self
                    .field
                    .exp(&domain_current_index, &T::from_u32(shift, self.field().p()));
                terms.push(self.field.mul(&q, &exp));
                power_map.insert(shift, exp);
            }

            for j in 0..usize::try_from(self.register_count).unwrap() {
                let bqv =
                    T::from_bytes_le(leaves[j].get(&current_index).unwrap(), self.field().p());
                terms.push(bqv.clone());
                let shift = transition_constraints_max_degree - boundary_quotient_degree_bounds[j];
                if let Some(exp) = power_map.get(&shift) {
                    terms.push(self.field.mul(&bqv, exp));
                } else {
                    let exp = self
                        .field
                        .exp(&domain_current_index, &T::from_u32(shift, self.field().p()));
                    terms.push(self.field.mul(&bqv, &exp));
                    power_map.insert(shift, exp);
                }
            }

            let mut combination = self.field.zero();
            for j in 0..terms.len() {
                combination = self
                    .field
                    .add(&combination, &self.field.mul(&terms[j], &weights[j]));
            }
            if combination != values[i] {
                panic!("invalid combination value");
            }
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use num_bigint::BigInt;

    use crate::{to_crypto_element, to_crypto_params, BigIntElement};

    use super::*;

    #[test]
    fn should_make_verify_stark_proof() {
        let p = to_crypto_params(BigIntElement(
            BigInt::from(1) + BigInt::from(407) * BigInt::from(2).pow(119),
        ));
        let g = to_crypto_element(
            BigIntElement(BigInt::from(85408008396924667383611388730472331217_u128)),
            &p,
        );
        let f = Rc::new(Field::new(g.clone()));

        let sequence_len = 40;
        let stark = Stark::new(&g.clone(), &f, 2, sequence_len, 32, 26, 2);

        let mut trace = Vec::new();
        trace.push(vec![f.bigint(2), f.bigint(3)]);
        trace.push(vec![f.bigint(4), f.bigint(9)]);
        while trace.len() < sequence_len.try_into().unwrap() {
            let e1 = &trace[trace.len() - 1][0];
            let e2 = &trace[trace.len() - 1][1];
            trace.push(vec![f.mul(&e1, &e1), f.mul(&e2, &e2)]);
        }

        let boundary_constraints = vec![
            (0, 0, f.bigint(2)),
            (0, 1, f.bigint(3)),
            (sequence_len - 1, 0, trace[trace.len() - 1][0].clone()),
            (sequence_len - 1, 1, trace[trace.len() - 1][1].clone()),
        ];

        let variables = MPolynomial::variables(1 + 2 * 2, &f);

        let _cycle_index = &variables[0];
        let prev_state = &variables[1..3];
        let next_state = &variables[3..];
        let mut transition_constraints = Vec::new();
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
        let ins = Instant::now();
        let proof = stark.prove(&trace, &transition_constraints, &boundary_constraints);
        println!("prove: {:?}", ins.elapsed());
        stark.verify(&proof, &transition_constraints, &boundary_constraints);
        println!("verify: {:?}", ins.elapsed());
    }
}
