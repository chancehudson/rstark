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

fn stark_(trace_len: u32, register_count: u32) -> Stark<BigIntElement> {
    let p = BigIntElement(BigInt::from(1) + BigInt::from(407) * BigInt::from(2).pow(119));
    let g = BigIntElement(BigInt::from(85408008396924667383611388730472331217_u128));
    let f = Rc::new(Field::new(p, g.clone()));

    Stark::<BigIntElement>::new(
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
    let input: ProveInput<BigIntElement> = serde_wasm_bindgen::from_value(input).unwrap();

    let register_count = input.trace[0].len();
    for i in 1..input.trace.len() {
        if input.trace[i].len() != register_count {
            log("inconsistent trace register count");
            panic!();
        }
    }
    let stark: Stark<BigIntElement> = stark_(
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
    let input: VerifyInput<BigIntElement> = serde_wasm_bindgen::from_value(input).unwrap();

    let stark = stark_(input.trace_len, input.register_count);

    let transition_constraints: Vec<MPolynomial<BigIntElement>> = input
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
    fn add(&self, v: &Self) -> Self;
    fn mul(&self, v: &Self) -> Self;
    fn sub(&self, v: &Self) -> Self;
    fn div(&self, v: &Self) -> Self;
    fn modd(&self, v: &Self) -> Self;
    fn neg(&self) -> Self;
    fn two() -> Self;
    fn one() -> Self;
    fn zero() -> Self;
    fn is_minus(&self) -> bool;
    fn extended_gcd(&self, v: &Self) -> Self;
    fn bits(&self) -> u64;
    fn from_u32(v: u32) -> Self;
    fn from_i32(v: i32) -> Self;
    fn to_u32(&self) -> u32;
    fn from_bytes_le(v: &[u8]) -> Self;
    fn to_bytes_le(&self) -> Vec<u8>;
    fn to_bytes_le_sized(&self) -> [u8; 32];
    fn modpow(&self, e: &Self, m: &Self) -> Self;
}

#[derive(Clone, PartialEq, Eq, Debug, PartialOrd, Hash, Deserialize)]
pub struct BigIntElement(pub BigInt);

impl FieldElement for BigIntElement {
    fn to_u32(&self) -> u32 {
        let (_, digits) = self.0.to_u32_digits();
        if digits.is_empty() {
            if digits.is_empty() {
                return 0;
            }
            panic!("invalid bigint digits len for u32 conversion");
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
    fn modpow(&self, e: &Self, m: &Self) -> Self {
        BigIntElement(self.0.modpow(&e.0, &m.0))
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
    fn add(&self, v: &Self) -> Self {
        BigIntElement(&self.0 + &v.0)
    }
    fn div(&self, v: &Self) -> Self {
        BigIntElement(&self.0 / &v.0)
    }
    fn mul(&self, v: &Self) -> Self {
        BigIntElement(&self.0 * &v.0)
    }

    fn sub(&self, v: &Self) -> Self {
        BigIntElement(&self.0 - &v.0)
    }

    fn modd(&self, v: &Self) -> Self {
        BigIntElement(&self.0 % &v.0)
    }

    fn neg(&self) -> Self {
        BigIntElement(&self.0 * BigInt::from(-1))
    }

    fn is_minus(&self) -> bool {
        self.0.sign() == Sign::Minus
    }
}
