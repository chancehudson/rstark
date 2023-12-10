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
use crypto_bigint::{Encoding, U256};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use serde::{Deserialize, Serialize};
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
    let p = to_crypto_params(BigIntElement(
        BigInt::from(1) + BigInt::from(407) * BigInt::from(2).pow(119),
    ));
    let g = to_crypto_element(
        BigIntElement(BigInt::from(85408008396924667383611388730472331217_u128)),
        &p,
    );
    let f = Rc::new(Field::new(g.clone()));

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

pub trait LegacyFieldElement: Eq + PartialEq + Clone + PartialOrd + Hash + Debug {
    fn add(&self, v: &Self, p: &Self) -> Self;
    fn mul(&self, v: &Self, p: &Self) -> Self;
    fn sub(&self, v: &Self, p: &Self) -> Self;
    fn div(&self, v: &Self, p: &Self) -> Self;
    fn modpow(&self, e: &Self, p: &Self) -> Self;
    fn modd(&self, p: &Self) -> Self;

    fn two() -> Self;
    fn one() -> Self;
    fn zero() -> Self;
    fn extended_gcd(&self, v: &Self) -> Self;
    fn bits(&self) -> u64;
    fn from_u32(v: u32) -> Self;
    fn from_i32(v: i32) -> Self;
    fn to_u32(&self) -> u32;
    fn from_bytes_le(v: &[u8]) -> Self;
    fn to_bytes_le(&self) -> Vec<u8>;
    fn to_bytes_le_sized(&self) -> [u8; 32];
}

#[derive(Clone, PartialEq, Eq, Debug, PartialOrd, Hash, Deserialize)]
pub struct BigIntElement(pub BigInt);

impl LegacyFieldElement for BigIntElement {
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
    fn from_i32(v: i32) -> Self {
        BigIntElement(BigInt::from(v))
    }
    fn from_u32(v: u32) -> Self {
        BigIntElement(BigInt::from(v))
    }
    fn bits(&self) -> u64 {
        self.0.bits()
    }
    fn extended_gcd(&self, v: &Self) -> Self {
        let res = self.0.extended_gcd(&v.0);
        if res.gcd != BigInt::from(1) {
            panic!("gcd != 1")
        } else {
            BigIntElement(res.x)
        }
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

pub trait FieldElement: Eq + PartialEq + Clone + PartialOrd + Hash + Debug {
    type ParamsType: Serialize + Debug;
    fn add(&self, v: &Self) -> Self;
    fn sub(&self, v: &Self) -> Self;
    fn mul(&self, v: &Self) -> Self;
    fn div(&self, v: &Self) -> Self;
    fn modpow(&self, e: &Self) -> Self;
    fn inv(&self) -> Self;
    fn neg(&self) -> Self;

    fn zero(p: &Self::ParamsType) -> Self;
    fn one(p: &Self::ParamsType) -> Self;
    fn two(p: &Self::ParamsType) -> Self;

    fn bits(&self) -> u64;
    fn from_u32(v: u32, p: &Self::ParamsType) -> Self;
    fn from_i32(v: i32, p: &Self::ParamsType) -> Self;
    fn to_u32(&self) -> u32;
    fn from_bytes_le(v: &[u8], p: &Self::ParamsType) -> Self;
    fn to_bytes_le(&self) -> Vec<u8>;
    fn to_bytes_le_sized(&self) -> [u8; 32];

    fn from_params(v: &Self::ParamsType) -> Self;
    fn get_params(&self) -> Self::ParamsType;
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct CryptoBigIntElement(pub DynResidue<4>);

impl<'de> Deserialize<'de> for CryptoBigIntElement {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::de::Deserializer<'de>,
    {
        let bytes = Vec::deserialize(deserializer)?;
        Ok(CryptoBigIntElement::from_bytes_le(
            &bytes[..],
            &ParamWrapper(DynResidueParams::new(&U256::ZERO)),
        ))
    }
}

impl PartialOrd for CryptoBigIntElement {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.0.retrieve().cmp(&other.0.retrieve()))
    }
}

impl Hash for CryptoBigIntElement {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.retrieve().hash(state)
    }
}

impl FieldElement for CryptoBigIntElement {
    type ParamsType = ParamWrapper;

    fn add(&self, v: &Self) -> Self {
        CryptoBigIntElement(self.0 + v.0)
    }
    fn sub(&self, v: &Self) -> Self {
        CryptoBigIntElement(self.0 - v.0)
    }
    fn mul(&self, v: &Self) -> Self {
        CryptoBigIntElement(self.0 * v.0)
    }
    fn div(&self, v: &Self) -> Self {
        CryptoBigIntElement(self.0 * v.0.invert().0)
    }
    fn modpow(&self, e: &Self) -> Self {
        CryptoBigIntElement(self.0.pow(&e.0.retrieve()))
    }
    fn inv(&self) -> Self {
        CryptoBigIntElement(self.0.invert().0)
    }
    fn neg(&self) -> Self {
        CryptoBigIntElement(self.0.neg())
    }

    fn zero(p: &ParamWrapper) -> Self {
        CryptoBigIntElement(DynResidue::new(&U256::ZERO, p.0))
    }
    fn one(p: &ParamWrapper) -> Self {
        CryptoBigIntElement(DynResidue::new(&U256::ONE, p.0))
    }
    fn two(p: &ParamWrapper) -> Self {
        CryptoBigIntElement(DynResidue::new(&U256::from_u32(2), p.0))
    }

    fn to_u32(&self) -> u32 {
        u32::from_le_bytes(self.to_bytes_le_sized()[..LIMBS].try_into().unwrap())
    }

    fn from_bytes_le(v: &[u8], p: &ParamWrapper) -> Self {
        CryptoBigIntElement(DynResidue::new(&U256::from_le_slice(v), p.0))
    }

    fn to_bytes_le(&self) -> Vec<u8> {
        self.0.retrieve().to_le_bytes().to_vec()
    }

    fn to_bytes_le_sized(&self) -> [u8; 32] {
        let mut extended = self.to_bytes_le();
        extended.resize(32, 0);
        extended.try_into().unwrap()
    }
    fn from_i32(v: i32, p: &ParamWrapper) -> Self {
        if v.is_negative() {
            let s = DynResidue::new(&U256::from_u32(-v as u32), p.0);
            // let s = -v
            // calculate (s/p + 1)*p - s
            CryptoBigIntElement(
                (s * Self::from_params(p).inv().0 + Self::one(p).0) * Self::from_params(p).0 - s,
            )
        } else {
            let val = DynResidue::new(&U256::from_u32(v as u32), p.0);
            CryptoBigIntElement(val)
        }
    }
    fn from_u32(v: u32, p: &ParamWrapper) -> Self {
        CryptoBigIntElement(DynResidue::new(&U256::from_u32(v), p.0))
    }
    fn bits(&self) -> u64 {
        self.0.retrieve().bits().try_into().unwrap()
    }
    fn from_params(v: &Self::ParamsType) -> Self {
        CryptoBigIntElement(DynResidue::new(v.0.modulus(), v.0))
    }
    fn get_params(&self) -> Self::ParamsType {
        ParamWrapper(*self.0.params())
    }
}

#[cfg(target_pointer_width = "16")]
const POINTER_WIDTH: usize = 16;

#[cfg(target_pointer_width = "32")]
const POINTER_WIDTH: usize = 32;

#[cfg(target_pointer_width = "64")]
const POINTER_WIDTH: usize = 64;

const LIMBS: usize = 256 / POINTER_WIDTH;
#[derive(Debug)]
pub struct ParamWrapper(pub DynResidueParams<LIMBS>);

impl Serialize for ParamWrapper {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.0.modulus().to_string())
    }
}

pub fn to_crypto_element(big_int: BigIntElement, p: &ParamWrapper) -> CryptoBigIntElement {
    let mut input = [0u8; U256::BYTES];
    let encoded = big_int.to_bytes_le();
    let l = encoded.len().min(U256::BYTES);
    input[..l].copy_from_slice(&encoded[..l]);
    CryptoBigIntElement(DynResidue::new(&U256::from_le_slice(&input), p.0))
}

pub fn to_crypto_params(big_int: BigIntElement) -> ParamWrapper {
    let mut input = [0u8; U256::BYTES];
    let encoded = big_int.to_bytes_le();
    let l = encoded.len().min(U256::BYTES);
    input[..l].copy_from_slice(&encoded[..l]);
    ParamWrapper(DynResidueParams::new(&U256::from_le_slice(&input)))
}
