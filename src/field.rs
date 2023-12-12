use crate::field_element::{FieldElement};
use rand::Rng;
use serde::Serialize;
use std::collections::HashMap;
use std::sync::RwLock;

#[derive(Serialize, Debug)]
pub struct Field<T: FieldElement> {
    p: T::ParamsType,
    g: T,
    // generator, size size keyed to contents
    group_cache: RwLock<HashMap<(T, u32), Vec<T>>>,
    // group size, offset keyed to contents
    coset_cache: RwLock<HashMap<(u32, T), Vec<T>>>,
    generator_cache: HashMap<u32, (T, T)>,
}

impl<T: FieldElement> Field<T> {
    pub fn new(g: T) -> Field<T> {
        let mut f = Field {
            p: g.get_params(),
            g,
            group_cache: RwLock::new(HashMap::new()),
            coset_cache: RwLock::new(HashMap::new()),
            generator_cache: HashMap::new(),
        };
        if T::from_params(&f.p).bits() <= 32 {
            // we're likely in the 101 field in tests
            // don't build the domain cache
            return f;
        }

        // build a cache of generators and inverted generators
        let mut start = 1;
        let mut sizes: Vec<u32> = Vec::new();
        let mut generators = Vec::new();
        for _ in 0..31 {
            start *= 2;
            generators.push(f.generator(T::from_u32(start, &f.p)));
            sizes.push(start);
        }
        let generators_inv = f.inv_batch(&generators);
        for i in 0..(generators.len()) {
            f.generator_cache
                .insert(sizes[i], (generators[i].clone(), generators_inv[i].clone()));
        }
        f
    }

    pub fn bigint(&self, val: i32) -> T {
        T::from_i32(val, self.p())
    }

    pub fn biguint(&self, val: u32) -> T {
        T::from_u32(val, self.p())
    }

    pub fn p(&self) -> &T::ParamsType {
        &self.p
    }

    pub fn g(&self) -> &T {
        &self.g
    }

    pub fn zero(&self) -> T {
        T::zero(self.p())
    }

    pub fn one(&self) -> T {
        T::one(self.p())
    }

    pub fn two(&self) -> T {
        T::two(self.p())
    }

    pub fn add(&self, v1: &T, v2: &T) -> T {
        v1.add(v2)
    }

    pub fn mul(&self, v1: &T, v2: &T) -> T {
        v1.mul(v2)
    }

    pub fn sub(&self, v1: &T, v2: &T) -> T {
        v1.sub(v2)
    }

    pub fn neg(&self, v: &T) -> T {
        v.neg()
    }

    pub fn div(&self, v1: &T, v2: &T) -> T {
        v1.div(v2)
    }

    // exponent should always be >= 0
    pub fn exp(&self, v: &T, e: &T) -> T {
        if e == &self.one() {
            v.clone()
        } else {
            v.modpow(e)
        }
    }

    pub fn generator_cache(&self, size: &u32) -> (T, T) {
        if let Some(v) = self.generator_cache.get(size) {
            return v.clone();
        }
        let g = self.generator(T::from_u32(*size, self.p()));
        let g_inv = self.inv(&g);
        (g, g_inv)
    }

    pub fn generator(&self, size: T) -> T {
        let transformed_p = T::from_params(self.p());
        let numer = &transformed_p.sub(&T::one(self.p()));
        let exp = numer.div(&size);
        if exp.mul(&size) != *numer {
            panic!("subgroup is not a divisor of field");
        }
        self.exp(&self.g, &exp)
    }

    pub fn inv(&self, v: &T) -> T {
        v.inv()
    }

    pub fn inv_batch(&self, values: &Vec<T>) -> Vec<T> {
        let mut last = T::one(self.p());
        let mut result = Vec::new();
        for v in values {
            result.push(last.clone());
            last = self.mul(&last, v);
        }
        last = self.inv(&last);
        for i in (0..values.len()).rev() {
            result[i] = self.mul(&result[i], &last);
            last = self.mul(&values[i], &last);
        }
        result
    }

    pub fn random(&self) -> T {
        let mut rng = rand::thread_rng();
        self.bigint(rng.gen())
    }

    pub fn sample(&self, input: T) -> T {
        // input.modd(self.p())
        input
    }

    pub fn sample_bytes(&self, input: &[u8]) -> T {
        let mut out = T::zero(self.p());
        for i in 0..16 {
            let i_b = T::from_u32(i, self.p());
            out = self.add(
                &out,
                &self.mul(
                    &T::from_u32(input[i as usize] as u32, self.p()),
                    &self.exp(&T::two(self.p()), &i_b),
                ),
            );
        }
        out
    }
    pub fn coset(&self, size: u32, offset: &T) -> Vec<T> {
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
        let generator = self.generator(T::from_u32(size, self.p()));
        // build, insert, and return
        let domain = self.domain(&generator, size);
        let mut coset = Vec::new();
        for v in domain {
            coset.push(self.mul(&v, offset));
        }
        let d = coset.clone();
        cache.insert((size, offset.clone()), coset);
        d
    }

    // retrieve the full domain, build if necessary
    pub fn domain(&self, generator: &T, size: u32) -> Vec<T> {
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
        let mut domain = Vec::new();
        domain.push(T::one(self.p()));
        for _ in 1..size {
            domain.push(self.mul(&domain[domain.len() - 1], generator));
        }
        let d = domain.clone();
        cache.insert((generator.clone(), size), domain);
        d
    }
}

#[cfg(test)]
mod tests {
    use crate::field_element::{ParamWrapper, CryptoBigIntElement, UC};
    use crypto_bigint::modular::runtime_mod::{DynResidueParams};

    use super::*;

    fn test_field() -> Field<CryptoBigIntElement> {
        let p = ParamWrapper(DynResidueParams::new(&UC::from_u128(101_u128)));
        let g = CryptoBigIntElement::from_u32(0, &p);
        Field::new(g)
    }

    #[test]
    fn should_make_bigint() {
        let f = test_field();
        assert_eq!(
            f.bigint(0),
            CryptoBigIntElement::from_u32(0, f.p())
        );
        assert_eq!(
            f.bigint(-29),
            CryptoBigIntElement::from_u32(72, f.p())
        );
        assert_eq!(
            f.bigint(-145),
            CryptoBigIntElement::from_u32(57, f.p())
        );
        assert_eq!(
            f.bigint(32),
            CryptoBigIntElement::from_u32(32, f.p())
        );
        assert_eq!(
            f.biguint(0),
            CryptoBigIntElement::from_u32(0, f.p())
        );
        assert_eq!(
            f.biguint(32),
            CryptoBigIntElement::from_u32(32, f.p())
        );
    }

    #[test]
    fn should_add_two_elements() {
        let f = test_field();

        let x = f.bigint(40);
        let y = f.bigint(90);

        assert_eq!(f.add(&x, &y), f.bigint(29));
    }

    #[test]
    fn should_add_neg_elements() {
        let f = test_field();

        let x = f.bigint(40);
        let y = f.bigint(90);

        assert_eq!(f.add(&x, &y), f.bigint(29));
    }

    #[test]
    fn should_mul_two_elements() {
        let f = test_field();

        let x = f.bigint(40);
        let y = f.bigint(90);

        assert_eq!(f.mul(&x, &y), f.bigint(65));
    }

    #[test]
    fn should_sub_two_elements() {
        let f = test_field();

        let x = f.bigint(2);
        let y = f.bigint(20);

        assert_eq!(f.sub(&x, &y), f.biguint(83));
    }

    #[test]
    fn should_get_generator() {
        let p = ParamWrapper(DynResidueParams::new(&UC::from_u128(3221225473_u128)));
        let g = CryptoBigIntElement::from_u32(5, &p);

        let f = Field::new(g);

        for i in 1..10 {
            let g = f.generator(f.biguint(u32::pow(2, i)));
            assert_eq!(f.exp(&g, &f.biguint(u32::pow(2, i))), f.bigint(1));
        }
    }

    #[test]
    fn should_get_inverse() {
        let p = ParamWrapper(DynResidueParams::new(&UC::from_u128(3221225473_u128)));
        let g = CryptoBigIntElement::from_u32(5, &p);
        let f = Field::new(g);

        for i in 1..99 {
            let v = f.biguint(i);
            let inv = f.inv(&v);
            assert_eq!(f.mul(&inv, &v), f.bigint(1));
        }
    }
}
