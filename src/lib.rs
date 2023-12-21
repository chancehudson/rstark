pub mod channel;
pub mod field;
pub mod fri;
pub mod mpolynomial;
pub mod polynomial;
pub mod stark;
pub mod tree;
pub mod field_element;

use crate::field::Field;
use crate::mpolynomial::MPolynomial;
use crate::stark::Stark;
use crate::field_element::{G, FieldElement, CryptoBigIntElement};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::rc::Rc;
use wasm_bindgen::prelude::*;

#[derive(Serialize, Deserialize)]
pub struct ProveInput<T: FieldElement> {
    trace: Vec<Vec<T>>,
    transition_constraints: Vec<HashMap<Vec<u32>, T>>,
    boundary: Vec<(u32, u32, T)>,
    perfect_zk: bool
}

#[derive(Serialize, Deserialize)]
pub struct VerifyInput<T: FieldElement> {
    trace_len: u32,
    register_count: u32,
    transition_constraints: Vec<HashMap<Vec<u32>, T>>,
    boundary: Vec<(u32, u32, T)>,
    perfect_zk: bool
}

fn stark_(trace_len: u32, register_count: u32, perfect_zk: bool) -> Stark<CryptoBigIntElement> {
    let f = Rc::new(Field::new(G));
    Stark::<CryptoBigIntElement>::new(
        &G,
        &f,
        register_count,
        trace_len,
        32, // expansion factor
        26, // colinearity tests
        2,  // constraint degree
        perfect_zk
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
        input.perfect_zk
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

    let stark = stark_(input.trace_len, input.register_count, input.perfect_zk);

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
