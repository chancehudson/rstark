use crate::channel::Channel;
use crate::field::Field;
use crate::polynomial::Polynomial;
use crate::tree::Tree;
use crate::FieldElement;
use std::collections::HashMap;
use std::rc::Rc;

pub struct FriOptions<T: FieldElement> {
    pub offset: T,
    pub omega: T,
    pub domain_len: u32,
    pub expansion_factor: u32,
    pub colinearity_test_count: u32,
}

pub struct Fri<T: FieldElement> {
    pub offset: T,
    pub omega: T,
    pub domain_len: u32,
    pub field: Rc<Field<T>>,
    pub expansion_factor: u32,
    pub colinearity_test_count: u32,
    domain: Vec<T>,
    round_count: u32,
}

impl<T: FieldElement> Fri<T> {
    pub fn new(options: &FriOptions<T>, field: &Rc<Field<T>>) -> Fri<T> {
        // calculate number of rounds
        let mut codeword_len = options.domain_len;
        let mut round_count = 0;
        while codeword_len > options.expansion_factor
            && 4 * options.colinearity_test_count < codeword_len
        {
            codeword_len /= 2;
            round_count += 1;
        }
        Fri {
            offset: options.offset.clone(),
            omega: options.omega.clone(),
            domain_len: options.domain_len,
            field: Rc::clone(field),
            expansion_factor: options.expansion_factor,
            colinearity_test_count: options.colinearity_test_count,
            domain: field.coset(options.domain_len, &options.offset),
            round_count,
        }
    }

    pub fn domain(&self) -> &Vec<T> {
        &self.domain
    }

    fn round_count(&self) -> u32 {
        self.round_count
    }

    pub fn prove(&self, codeword: &Vec<T>, channel: &mut Channel) -> Vec<u32> {
        if self.domain_len != u32::try_from(codeword.len()).unwrap() {
            panic!("initial codeword does not match domain len");
        }
        let codewords = self.commit(codeword, channel);
        let top_indices = self.sample_indices(
            &channel.prover_hash(),
            codewords[1].len().try_into().unwrap(),
            codewords[codewords.len() - 1].len().try_into().unwrap(),
            self.colinearity_test_count,
        );
        let mut indices: Vec<u32> = top_indices.clone();
        let codeword_trees: Vec<Tree<T>> = codewords
            .iter()
            .map(|word| {
                Tree::build(
                    &word
                        .iter()
                        .map(|t| t.to_bytes_le_sized())
                        .collect::<Vec<[u8; 32]>>(),
                )
            })
            .collect();
        for i in 0..(codewords.len() - 1) {
            indices = indices
                .iter()
                .map(|index| index % u32::try_from(codewords[i].len() >> 1).unwrap())
                .collect();
            self.query(
                &codewords[i],
                &codewords[i + 1],
                &indices,
                channel,
                &codeword_trees[i],
                &codeword_trees[i + 1],
            );
        }
        top_indices
    }

    fn query(
        &self,
        current_codeword: &Vec<T>,
        next_codeword: &[T],
        indices_c: &[u32],
        channel: &mut Channel,
        current_codeword_tree: &Tree<T>,
        next_codeword_tree: &Tree<T>,
    ) {
        let indices_a: Vec<u32> = indices_c.to_vec();
        let indices_b: Vec<u32> = indices_c
            .iter()
            .map(|val| val + ((current_codeword.len() >> 1) as u32))
            .collect();
        for i in 0..usize::try_from(self.colinearity_test_count).unwrap() {
            channel.push(&[
                current_codeword[usize::try_from(indices_a[i]).unwrap()].to_bytes_le_sized(),
                current_codeword[usize::try_from(indices_b[i]).unwrap()].to_bytes_le_sized(),
                next_codeword[usize::try_from(indices_c[i]).unwrap()].to_bytes_le_sized(),
            ]);
        }
        for i in 0..usize::try_from(self.colinearity_test_count).unwrap() {
            channel.push(&current_codeword_tree.open(indices_a[i]).0);
            channel.push(&current_codeword_tree.open(indices_b[i]).0);
            channel.push(&next_codeword_tree.open(indices_c[i]).0);
        }
    }

    fn commit(&self, codeword: &[T], channel: &mut Channel) -> Vec<Vec<T>> {
        let mut codewords = Vec::new();
        let mut codeword = codeword.to_owned();
        let two_inv = self.field.inv(&self.field.two());

        // invert the entire domain using repeated multiplications
        // e.g. 1/4 = (1/2) * (1/2)
        // 1/x^2 = (1/x) * (1/x)
        let inv_offset = self.field.inv(&self.offset);
        let inv_offset_domain = self
            .field
            .domain(&inv_offset, 2_u32.pow(self.round_count()));

        let inv_omega = self.field.inv(&self.omega);
        let inv_domain = self.field.domain(&inv_omega, self.domain_len);

        let mut exp: usize = 1;

        for x in 0..self.round_count() {
            let root = Tree::commit_elements(&codeword);
            channel.push_single(&root);
            if x == self.round_count() - 1 {
                break;
            }
            codewords.push(codeword.clone());
            // now split the last codeword and fold into a set
            // of points from a polynomial of half the degree
            // of the previous codewords, similar to a FFT
            let alpha = self.field.sample_bytes(&channel.prover_hash());

            let next_len = codeword.len() >> 1;
            codeword = codeword[0..next_len]
                .iter()
                .enumerate()
                .map(|(index, val)| {
                    let inv_omega = self
                        .field
                        .mul(&inv_offset_domain[exp], &inv_domain[exp * index]);
                    // ( (one + alpha / (offset * (omega^i)) ) * codeword[i]
                    let a = self.field.mul(
                        val,
                        &self
                            .field
                            .add(&self.field.one(), &self.field.mul(&alpha, &inv_omega)),
                    );
                    //  (one - alpha / (offset * (omega^i)) ) * codeword[len(codeword)//2 + i] ) for i in range(len(codeword)//2)]
                    let b = self.field.mul(
                        &self
                            .field
                            .sub(&self.field.one(), &self.field.mul(&alpha, &inv_omega)),
                        &codeword[(codeword.len() >> 1) + index],
                    );
                    self.field.mul(&two_inv, &self.field.add(&a, &b))
                })
                .collect();

            exp *= 2;
        }

        channel.push(
            &codeword
                .iter()
                .map(|t| t.to_bytes_le_sized())
                .collect::<Vec<[u8; 32]>>(),
        );
        codewords.push(codeword);
        codewords
    }

    fn sample_indices(
        &self,
        seed: &[u8; 32],
        size: u32,
        reduced_size: u32,
        count: u32,
    ) -> Vec<u32> {
        if count > 2 * reduced_size {
            panic!("not enough entropy");
        }
        if count > reduced_size {
            panic!("cannot sample more indices than available");
        }

        let mut indices: Vec<u32> = Vec::new();
        let mut reduced_indices: HashMap<u32, bool> = HashMap::new();
        let mut counter: u32 = 0;
        while indices.len() < (count as usize) {
            let mut hasher = blake3::Hasher::new();
            hasher.update(seed);
            hasher.update(&T::from_u32(counter, self.field.p()).to_bytes_le());
            let v = T::from_bytes_le(hasher.finalize().as_bytes(), self.field.p()).to_u32();
            let index = v % size;
            let reduced_index = index % reduced_size;
            counter += 1;
            reduced_indices.entry(reduced_index).or_insert_with(|| {
                indices.push(index);
                true
            });
        }
        indices
    }

    pub fn verify(&self, channel: &mut Channel) -> Vec<(u32, T)> {
        let mut out = Vec::new();
        let mut offset = self.offset.clone();

        let mut roots = Vec::new();
        let mut alphas = Vec::new();

        for _ in 0..self.round_count() {
            roots.push(channel.pull_root());
            alphas.push(self.field.sample_bytes(&channel.verifier_hash()));
        }

        let last_codeword = channel.pull_path();
        if roots[roots.len() - 1] != Tree::<T>::commit(&last_codeword) {
            panic!("last codeword root mismatch");
        }

        let degree: usize =
            (last_codeword.len() / usize::try_from(self.expansion_factor).unwrap()) - 1;
        let mut last_offset = offset.clone();
        for _ in 0..(self.round_count() - 1) {
            last_offset = self.field.mul(&last_offset, &last_offset);
        }

        let omega_domain = self.field.domain(&self.omega, self.domain_len);
        let omega_start_index = 2_usize.pow(self.round_count() - 1);
        if self.field.inv(&omega_domain[omega_start_index])
            != self.field.exp(
                &omega_domain[omega_start_index],
                &T::from_u32((last_codeword.len() - 1) as u32, self.field.p()),
            )
        {
            panic!("omega order incorrect");
        }

        let last_domain = last_codeword
            .iter()
            .enumerate()
            .map(|(index, _)| {
                self.field
                    .mul(&last_offset, &omega_domain[omega_start_index * index])
            })
            .collect();

        let poly = Polynomial::interpolate_fft(
            &last_domain,
            &last_codeword
                .iter()
                .map(|v| T::from_bytes_le(v, self.field.p()))
                .collect(),
            &self.field,
        );
        for i in 0..last_domain.len() {
            if poly.eval(&last_domain[i]) != T::from_bytes_le(&last_codeword[i], self.field.p()) {
                panic!("interpolated polynomial is incorrect");
            }
        }
        if poly.degree() > degree {
            panic!("last codeword does not match polynomial of low enough degree");
        }

        let top_indices = self.sample_indices(
            &channel.verifier_hash(),
            self.domain_len >> 1,
            self.domain_len >> (self.round_count() - 1),
            self.colinearity_test_count,
        );
        let mut colinearity_x_vals = Vec::new();
        let mut colinearity_y_vals = Vec::new();
        let mut exp = 1;
        for i in 0..usize::try_from(self.round_count() - 1).unwrap() {
            let indices_c: Vec<u32> = top_indices
                .iter()
                .map(|val| val % (self.domain_len >> (i + 1)))
                .collect();
            let indices_a: Vec<u32> = indices_c.clone();
            let indices_b: Vec<u32> = indices_a
                .iter()
                .map(|val| val + (self.domain_len >> (i + 1)))
                .collect();

            let mut aa = Vec::new();
            let mut bb = Vec::new();
            let mut cc = Vec::new();
            for j in 0..usize::try_from(self.colinearity_test_count).unwrap() {
                let y_points_msg = channel.pull_path();
                let ay = T::from_bytes_le(&y_points_msg[0], self.field.p());
                let by = T::from_bytes_le(&y_points_msg[1], self.field.p());
                let cy = T::from_bytes_le(&y_points_msg[2], self.field.p());
                aa.push(ay.clone());
                bb.push(by.clone());
                cc.push(cy.clone());
                if i == 0 {
                    out.push((indices_a[j], ay.clone()));
                    out.push((indices_b[j], by.clone()));
                }

                let index_a_usize = usize::try_from(indices_a[j]).unwrap();
                let index_b_usize = usize::try_from(indices_b[j]).unwrap();
                let ax = self
                    .field
                    .mul(&offset, &omega_domain[exp * index_a_usize])
                    .clone();
                let bx = self
                    .field
                    .mul(&offset, &omega_domain[exp * index_b_usize])
                    .clone();

                let cx = alphas[i].clone();

                colinearity_x_vals.push(vec![ax, bx, cx]);
                colinearity_y_vals.push(vec![ay, by, cy]);
            }

            if !Polynomial::test_colinearity_batch(
                &colinearity_x_vals,
                &colinearity_y_vals,
                &self.field,
            ) {
                panic!("colinearity test failed");
            }

            for j in 0..usize::try_from(self.colinearity_test_count).unwrap() {
                Tree::<T>::verify(
                    &roots[i],
                    indices_a[j],
                    &channel.pull_path(),
                    &aa[j].to_bytes_le_sized(),
                );
                Tree::<T>::verify(
                    &roots[i],
                    indices_b[j],
                    &channel.pull_path(),
                    &bb[j].to_bytes_le_sized(),
                );
                Tree::<T>::verify(
                    &roots[i + 1],
                    indices_c[j],
                    &channel.pull_path(),
                    &cc[j].to_bytes_le_sized(),
                );
            }

            exp *= 2;
            offset = self.field.mul(&offset, &offset);
        }

        out
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;

    use crate::{to_crypto_element, to_crypto_params, BigIntElement};

    use super::*;

    #[test]
    fn should_make_verify_fri_proof() {
        let mut channel = Channel::new();
        let p = to_crypto_params(BigIntElement(
            BigInt::from(1) + BigInt::from(407) * BigInt::from(2).pow(119),
        ));
        let g = to_crypto_element(
            BigIntElement(BigInt::from(85408008396924667383611388730472331217_u128)),
            &p,
        );
        let f = Rc::new(Field::new(g.clone()));

        let domain_size: u32 = 8192;
        let domain_g = f.generator(f.biguint(domain_size));

        let fri = Fri::new(
            &FriOptions {
                offset: g.clone(),
                omega: domain_g.clone(),
                domain_len: domain_size,
                expansion_factor: 2,
                colinearity_test_count: 10,
            },
            &f,
        );

        let mut poly = Polynomial::new(&f);
        poly.term(&f.bigint(3), 2);
        let mut points = Vec::new();
        for i in fri.domain() {
            points.push(poly.eval(i));
        }
        fri.prove(&points, &mut channel);
        fri.verify(&mut channel);
    }
}
