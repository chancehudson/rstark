use crate::FieldElement;
use rand::Rng;
use serde::Serialize;
use std::collections::HashMap;
use std::sync::RwLock;

#[derive(Serialize, Debug)]
pub struct Field<T: FieldElement> {
    p: T,
    g: T,
    // generator, size size keyed to contents
    group_cache: RwLock<HashMap<(T, u32), Vec<T>>>,
    // group size, offset keyed to contents
    coset_cache: RwLock<HashMap<(u32, T), Vec<T>>>,
    generator_cache: HashMap<u32, (T, T)>,
}

impl<T: FieldElement> Field<T> {
    pub fn new(p: T, g: T) -> Field<T> {
        let mut f = Field {
            p,
            g,
            group_cache: RwLock::new(HashMap::new()),
            coset_cache: RwLock::new(HashMap::new()),
            generator_cache: HashMap::new(),
        };
        if f.p.bits() <= 32 {
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
            generators.push(f.generator(T::from_u32(start)));
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
        self.modd(T::from_i32(val))
    }

    pub fn biguint(&self, val: u32) -> T {
        self.modd(T::from_u32(val))
    }

    pub fn p(&self) -> &T {
        &self.p
    }

    pub fn g(&self) -> &T {
        &self.g
    }

    pub fn zero() -> T {
        T::zero()
    }

    pub fn one() -> T {
        T::one()
    }

    pub fn modd(&self, v: T) -> T {
        if v.is_minus() {
            // v + p*((-v/p)+1)
            return v.add(&self.p.mul(&v.neg().div(&self.p).add(&T::one())));
        }
        v.modd(&self.p)
    }

    pub fn add(&self, v1: &T, v2: &T) -> T {
        self.modd(v1.add(v2))
    }

    pub fn ladd(&self, v1: &T, v2: &T) -> T {
        v1.add(v2)
    }

    pub fn mul(&self, v1: &T, v2: &T) -> T {
        self.modd(v1.mul(v2))
    }

    pub fn lmul(&self, v1: &T, v2: &T) -> T {
        v1.mul(v2)
    }

    pub fn sub(&self, v1: &T, v2: &T) -> T {
        self.modd(v1.sub(v2))
    }

    pub fn lsub(&self, v1: &T, v2: &T) -> T {
        v1.sub(v2)
    }

    pub fn neg(&self, v: &T) -> T {
        self.modd(self.p.sub(v))
    }

    pub fn div(&self, v1: &T, v2: &T) -> T {
        self.mul(v1, &self.inv(v2))
    }

    // exponent should always be >= 0
    pub fn exp(&self, v: &T, e: &T) -> T {
        if e == &Field::one() {
            v.clone()
        } else {
            self.modd(v.clone()).modpow(e, &self.p)
        }
    }

    pub fn generator_cache(&self, size: &u32) -> (T, T) {
        if let Some(v) = self.generator_cache.get(size) {
            return v.clone();
        }
        let g = self.generator(T::from_u32(*size));
        let g_inv = self.inv(&g);
        (g, g_inv)
    }

    pub fn generator(&self, size: T) -> T {
        if size >= self.p {
            panic!("requested subgroup is larger than field");
        }
        let numer = &self.p.sub(&T::one());
        let exp = numer.div(&size);
        if exp.mul(&size) != *numer {
            panic!("subgroup is not a divisor of field");
        }
        self.exp(&self.g, &exp)
    }

    pub fn inv(&self, v: &T) -> T {
        let e_gcd = v.clone().extended_gcd(&self.p);
        self.add(&e_gcd, &self.p)
    }

    pub fn inv_batch(&self, values: &Vec<T>) -> Vec<T> {
        let mut last = T::one();
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
        self.modd(input)
    }

    pub fn sample_bytes(&self, input: &[u8]) -> T {
        let mut out = T::zero();
        for i in 0..16 {
            let i_b = T::from_u32(i);
            out = self.add(
                &out,
                &self.mul(
                    &T::from_u32(input[i as usize] as u32),
                    &self.exp(&T::two(), &i_b),
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
        let generator = self.generator(T::from_u32(size));
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
        domain.push(T::one());
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
    use num_bigint::{BigInt, Sign};

    use crate::BigIntElement;

    use super::*;

    #[test]
    fn should_make_bigint() {
        let p = BigIntElement(BigInt::from(101));
        let f = Field::new(p, BigIntElement(BigInt::from(0)));
        assert_eq!(
            f.bigint(0),
            BigIntElement(BigInt::new(Sign::Minus, vec!(0)))
        );
        assert_eq!(
            f.bigint(-29),
            BigIntElement(BigInt::new(Sign::Plus, vec!(72)))
        );
        assert_eq!(
            f.bigint(32),
            BigIntElement(BigInt::new(Sign::Plus, vec!(32)))
        );
        assert_eq!(
            f.biguint(0),
            BigIntElement(BigInt::new(Sign::Minus, vec!(0)))
        );
        assert_eq!(
            f.biguint(32),
            BigIntElement(BigInt::new(Sign::Plus, vec!(32)))
        );
    }

    #[test]
    fn should_add_two_elements() {
        let p = BigIntElement(BigInt::from(101));
        let g = BigIntElement(BigInt::from(0));
        let f = Field::new(p, g);

        let x = f.bigint(40);
        let y = f.bigint(90);

        assert_eq!(f.add(&x, &y), f.bigint(29));
    }

    #[test]
    fn should_add_neg_elements() {
        let p = BigIntElement(BigInt::from(101));
        let g = BigIntElement(BigInt::from(0));
        let f = Field::new(p, g);

        let x = f.bigint(40);
        let y = f.bigint(90);

        assert_eq!(f.add(&x, &y), f.bigint(29));
    }
    #[test]
    fn should_mul_two_elements() {
        let p = BigIntElement(BigInt::from(101));
        let g = BigIntElement(BigInt::from(0));
        let f = Field::new(p, g);

        let x = f.bigint(40);
        let y = f.bigint(90);

        assert_eq!(f.mul(&x, &y), f.bigint(65));
    }

    #[test]
    fn should_sub_two_elements() {
        let p = BigIntElement(BigInt::from(101));
        let g = BigIntElement(BigInt::from(0));
        let f = Field::new(p, g);

        let x = f.bigint(2);
        let y = f.bigint(20);

        assert_eq!(f.sub(&x, &y), f.biguint(83));
    }

    #[test]
    fn should_get_generator() {
        let p = BigIntElement(BigInt::from(3221225473_u32));
        let f_g = BigIntElement(BigInt::from(5));
        let f = Field::new(p, f_g);

        for i in 1..10 {
            let g = f.generator(f.biguint(u32::pow(2, i)));
            assert_eq!(f.exp(&g, &f.biguint(u32::pow(2, i))), f.bigint(1));
        }
    }

    #[test]
    fn should_get_inverse() {
        let p = BigIntElement(BigInt::from(3221225473_u32));
        let f_g = BigIntElement(BigInt::from(5));
        let f = Field::new(p, f_g);

        for i in 1..99 {
            let v = f.biguint(i);
            let inv = f.inv(&v);
            assert_eq!(f.mul(&inv, &v), Field::one());
        }
    }
}
