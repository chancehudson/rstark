pub mod channel;
pub mod field;
pub mod fri;
pub mod mpolynomial;
pub mod polynomial;
pub mod stark;
pub mod tree;

use crate::field::Field;
use crate::mpolynomial::MPolynomial;
use crate::stark::Stark;
use crypto_bigint::modular::runtime_mod::{DynResidue, DynResidueParams};
use crypto_bigint::{Encoding, NonZero, U256};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use serde::{Deserialize, Deserializer, Serialize};
use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;
use std::rc::Rc;
use wasm_bindgen::prelude::*;
#[derive(Serialize, Deserialize)]
pub struct ProveInput<T: FieldElement> {
    trace: Vec<Vec<T>>,
    transition_constraints: Vec<HashMap<Vec<u32>, T>>,
    boundary: Vec<(u32, u32, T)>,
}

#[derive(Serialize, Deserialize)]
pub struct VerifyInput<T: FieldElement> {
    trace_len: u32,
    register_count: u32,
    transition_constraints: Vec<HashMap<Vec<u32>, T>>,
    boundary: Vec<(u32, u32, T)>,
}

fn stark_(trace_len: u32, register_count: u32) -> Stark<CryptoBigIntElement> {
    let p = CryptoBigIntElement(to_u256(
        BigInt::from(1) + BigInt::from(407) * BigInt::from(2).pow(119),
    ));
    let g = CryptoBigIntElement(to_u256(BigInt::from(
        85408008396924667383611388730472331217_u128,
    )));
    let f = Rc::new(Field::new(p, g.clone()));

    Stark::<CryptoBigIntElement>::new(
        &g,
        &f,
        register_count,
        trace_len,
        32, // expansion factor
        26, // colinearity tests
        2,  // constraint degree
    )
}

#[wasm_bindgen]
pub fn prove(input: JsValue) -> String {
    let input: ProveInput<CryptoBigIntElement> = serde_wasm_bindgen::from_value(input).unwrap();

    let register_count = input.trace[0].len();
    for i in 1..input.trace.len() {
        if input.trace[i].len() != register_count {
            log("inconsistent trace register count");
            panic!();
        }
    }
    let stark: Stark<CryptoBigIntElement> = stark_(
        input.trace.len().try_into().unwrap(),
        register_count.try_into().unwrap(),
    );

    let transition_constraints = input
        .transition_constraints
        .iter()
        .map(|v| MPolynomial::from_map(v, stark.field()))
        .collect();
    let boundary_constraints = input
        .boundary
        .iter()
        .map(|(v1, v2, v3)| (*v1, *v2, v3.clone()))
        .collect();
    let trace = input.trace.iter().map(|v| v.to_vec()).collect();

    stark.prove(&trace, &transition_constraints, &boundary_constraints)
}

#[wasm_bindgen]
pub fn verify(proof: String, input: JsValue) {
    let input: VerifyInput<CryptoBigIntElement> = serde_wasm_bindgen::from_value(input).unwrap();

    let stark = stark_(input.trace_len, input.register_count);

    let transition_constraints: Vec<MPolynomial<CryptoBigIntElement>> = input
        .transition_constraints
        .iter()
        .map(|v| MPolynomial::from_map(v, stark.field()))
        .collect();
    let boundary_constraints = input
        .boundary
        .iter()
        .map(|(v1, v2, v3)| (*v1, *v2, v3.clone()))
        .collect();
    stark.verify(&proof, &transition_constraints, &boundary_constraints);
}

#[wasm_bindgen]
extern "C" {
    // Use `js_namespace` here to bind `console.log(..)` instead of just
    // `log(..)`
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);

    // The `console.log` is quite polymorphic, so we can bind it with multiple
    // signatures. Note that we need to use `js_name` to ensure we always call
    // `log` in JS.
    #[wasm_bindgen(js_namespace = console, js_name = log)]
    fn log_u32(a: u32);

    // Multiple arguments too!
    #[wasm_bindgen(js_namespace = console, js_name = log)]
    fn log_many(a: &str, b: &str);
}

pub trait FieldElement: Eq + PartialEq + Clone + PartialOrd + Hash + Debug {
    fn add(&self, v: &Self, p: &Self) -> Self;
    fn mul(&self, v: &Self, p: &Self) -> Self;
    fn sub(&self, v: &Self, p: &Self) -> Self;
    fn div(&self, v: &Self, p: &Self) -> Self;
    fn modpow(&self, e: &Self, p: &Self) -> Self;
    fn modd(&self, p: &Self) -> Self;

    fn two() -> Self;
    fn one() -> Self;
    fn zero() -> Self;
    fn inv(&self, p: &Self) -> Self;
    fn bits(&self) -> u64;
    fn from_u32(v: u32) -> Self;
    fn from_i32(v: i32, p: &Self) -> Self;
    fn to_u32(&self) -> u32;
    fn from_bytes_le(v: &[u8]) -> Self;
    fn to_bytes_le(&self) -> Vec<u8>;
    fn to_bytes_le_sized(&self) -> [u8; 32];
}

#[derive(Clone, PartialEq, Eq, Debug, PartialOrd, Hash, Deserialize)]
pub struct BigIntElement(pub BigInt);

impl FieldElement for BigIntElement {
    fn to_u32(&self) -> u32 {
        let (_, digits) = self.0.to_u32_digits();
        if digits.is_empty() {
            return 0;
        }
        digits[0]
    }

    fn from_bytes_le(bytes: &[u8]) -> Self {
        BigIntElement(BigInt::from_bytes_le(Sign::Plus, bytes))
    }
    fn to_bytes_le(&self) -> Vec<u8> {
        self.0.to_bytes_le().1
    }
    fn to_bytes_le_sized(&self) -> [u8; 32] {
        let mut extended = self.to_bytes_le();
        extended.resize(32, 0);
        extended.try_into().unwrap()
    }
    fn from_i32(v: i32, p: &Self) -> Self {
        assert!(p.0.sign() == Sign::Plus);
        if v.is_negative() {
            BigIntElement((-v / &p.0 + 1) * &p.0 + v)
        } else {
            BigIntElement(v % &p.0)
        }
    }
    fn from_u32(v: u32) -> Self {
        BigIntElement(BigInt::from(v))
    }
    fn bits(&self) -> u64 {
        self.0.bits()
    }
    fn inv(&self, p: &Self) -> Self {
        let res = self.0.extended_gcd(&p.0);
        if res.gcd != BigInt::from(1) {
            panic!("gcd != 1");
        }
        BigIntElement(res.x).add(p, p)
    }
    fn two() -> Self {
        BigIntElement(BigInt::from(2))
    }
    fn one() -> Self {
        BigIntElement(BigInt::from(1))
    }
    fn zero() -> Self {
        BigIntElement(BigInt::from(0))
    }
    fn add(&self, v: &Self, p: &Self) -> Self {
        BigIntElement(&self.0 + &v.0).modd(p)
    }
    fn div(&self, v: &Self, p: &Self) -> Self {
        BigIntElement(&self.0 / &v.0).modd(p)
    }
    fn mul(&self, v: &Self, p: &Self) -> Self {
        BigIntElement(&self.0 * &v.0).modd(p)
    }
    fn sub(&self, v: &Self, p: &Self) -> Self {
        BigIntElement(&self.0 - &v.0).modd(p)
    }
    fn modd(&self, p: &Self) -> Self {
        assert!(p.0.sign() == Sign::Plus);
        if self.0.sign() == Sign::Minus {
            BigIntElement(((-&self.0) / &p.0 + 1) * &p.0 + &self.0)
        } else {
            BigIntElement(&self.0 % &p.0)
        }
    }
    fn modpow(&self, e: &Self, p: &Self) -> Self {
        BigIntElement(self.0.modpow(&e.0, &p.0))
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Hash)]
pub struct CryptoBigIntElement(pub U256);

impl Debug for CryptoBigIntElement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("check")
            .field("x", &to_biguint(&self.0).to_string())
            .finish()
    }
}
impl<'de> Deserialize<'de> for CryptoBigIntElement {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        Ok(CryptoBigIntElement(U256::from_le_hex(&s)))
    }
}

impl FieldElement for CryptoBigIntElement {
    fn add(&self, v: &Self, p: &Self) -> Self {
        CryptoBigIntElement(self.0.wrapping_add(&v.0).rem(&NonZero::new(p.0).unwrap()))
    }

    fn mul(&self, v: &Self, p: &Self) -> Self {
        CryptoBigIntElement(self.0.wrapping_mul(&v.0).rem(&NonZero::new(p.0).unwrap()))
    }

    fn sub(&self, v: &Self, p: &Self) -> Self {
        CryptoBigIntElement(self.0.sub_mod(&v.0, &p.0))
    }

    fn div(&self, v: &Self, p: &Self) -> Self {
        let s = self.0.checked_div(&v.0).unwrap();
        CryptoBigIntElement(s.rem(&NonZero::new(p.0).unwrap()))
    }

    fn modpow(&self, e: &Self, p: &Self) -> Self {
        // (self ^ exponent) mod modulus
        let params = DynResidueParams::new(&p.0);
        let a_m = DynResidue::new(&self.0, params);
        CryptoBigIntElement(a_m.pow(&e.0).retrieve())
    }

    fn modd(&self, p: &Self) -> Self {
        // assert!(p.0.sign() == Sign::Plus);
        // if self.0.sign() == Sign::Minus {
        //     BigIntElement(((-&self.0) / &p.0 + 1) * &p.0 + &self.0)
        // } else {
        //     BigIntElement(&self.0 % &p.0)
        // }

        CryptoBigIntElement(self.0.rem(&NonZero::new(p.0).unwrap()))
    }

    fn two() -> Self {
        CryptoBigIntElement(U256::from_u32(2))
    }

    fn one() -> Self {
        CryptoBigIntElement(U256::ONE)
    }

    fn zero() -> Self {
        CryptoBigIntElement(U256::ZERO)
    }

    fn inv(&self, p: &Self) -> Self {
        CryptoBigIntElement(self.0.inv_mod(&p.0).0)
    }
    fn bits(&self) -> u64 {
        self.0.bits().try_into().unwrap()
    }

    fn from_u32(v: u32) -> Self {
        CryptoBigIntElement(U256::from_u32(v))
    }

    fn from_i32(v: i32, p: &Self) -> Self {
        if v.is_negative() {
            let s = U256::from_u32(-v as u32);
            CryptoBigIntElement(
                s.wrapping_div(&p.0)
                    .wrapping_add(&U256::ONE)
                    .wrapping_mul(&p.0.wrapping_sub(&s)),
            )
        } else {
            CryptoBigIntElement(U256::from_u32(v as u32).rem(&NonZero::new(p.0).unwrap()))
        }
    }

    fn to_u32(&self) -> u32 {
        u32::from_le_bytes(self.to_bytes_le_sized()[..4].try_into().unwrap())
    }

    fn from_bytes_le(v: &[u8]) -> Self {
        CryptoBigIntElement(U256::from_le_slice(v))
    }

    fn to_bytes_le(&self) -> Vec<u8> {
        self.0.to_le_bytes().to_vec()
    }

    fn to_bytes_le_sized(&self) -> [u8; 32] {
        let mut extended = self.to_bytes_le();
        extended.resize(32, 0);
        extended.try_into().unwrap()
    }
}

pub fn to_u256(big_int: BigInt) -> U256 {
    let mut input = [0u8; U256::BYTES];
    let encoded = big_int.to_bytes_le();
    let l = encoded.1.len().min(U256::BYTES);
    input[..l].copy_from_slice(&encoded.1[..l]);

    U256::from_le_slice(&input)
}

fn to_biguint(uint: &U256) -> BigInt {
    BigInt::from_bytes_le(Sign::Plus, uint.to_le_bytes().as_ref())
}
