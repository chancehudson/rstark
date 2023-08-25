pub mod channel;
pub mod field;
pub mod polynomial;
pub mod tree;
pub mod mpolynomial;
pub mod fri;
pub mod stark;

use crate::field::Field;
use crate::stark::Stark;
use crate::mpolynomial::MPolynomial;
use std::rc::Rc;
use num_bigint::{BigInt, BigUint, ToBigInt};
use serde::{Serialize, Deserialize};
use std::collections::HashMap;

use wasm_bindgen::prelude::*;

#[derive(Serialize, Deserialize)]
pub struct ProveInput {
  trace: Vec<Vec<u128>>,
  transition_constraints: Vec<HashMap<Vec<u32>, u128>>,
  boundary: Vec<(u32, u32, u128)>
}

#[derive(Serialize, Deserialize)]
pub struct VerifyInput {
  trace_len: u32,
  register_count: u32,
  transition_constraints: Vec<HashMap<Vec<u32>, u128>>,
  boundary: Vec<(u32, u32, u128)>
}

fn stark_(trace_len: u32, register_count: u32) -> Stark {
  let p = 1 + 407 * 2_u128.pow(119);
  let g = 85408008396924667383611388730472331217_u128;
  let f = Rc::new(Field::new(p, g.clone()));

  Stark::new(
    &g.clone(),
    &f,
    register_count,
    trace_len,
    32, // expansion factor
    26, // colinearity tests
    2 // constraint degree
  )
}

#[wasm_bindgen]
pub fn prove(input: JsValue) -> String {
  let input: ProveInput = serde_wasm_bindgen::from_value(input).unwrap();

  let register_count = input.trace[0].len();
  for i in 1..input.trace.len() {
    if input.trace[i].len() != register_count {
      log("inconsistent trace register count");
      panic!();
    }
  }
  let stark = stark_(input.trace.len().try_into().unwrap(), register_count.try_into().unwrap());

  let transition_constraints = input.transition_constraints
    .iter()
    .map(|v| MPolynomial::from_map(&v, stark.field()))
    .collect();
  let boundary_constraints = input.boundary;
    // .iter()
    // .map(|(v1, v2, v3)| (v1.clone(), v2.clone(), v3.clone()))
    // .collect();
  let trace = input.trace;

  stark.prove(&trace, &transition_constraints, &boundary_constraints)
}

#[wasm_bindgen]
pub fn verify(proof: String, input: JsValue) {
  let input: VerifyInput = serde_wasm_bindgen::from_value(input).unwrap();

  let stark = stark_(input.trace_len, input.register_count);

  let transition_constraints: Vec<MPolynomial> = input.transition_constraints
    .iter()
    .map(|v| MPolynomial::from_map(&v, stark.field()))
    .collect();
  let boundary_constraints: Vec<(u32, u32, u128)> = input.boundary;
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
