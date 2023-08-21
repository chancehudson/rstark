mod channel;
mod field;
mod polynomial;
mod tree;
mod fri;
mod mpolynomial;
mod stark;

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
  trace: Vec<Vec<BigUint>>,
  transition_constraints: Vec<HashMap<Vec<u32>, BigUint>>,
  boundary: Vec<(u32, u32, BigUint)>
}

#[derive(Serialize, Deserialize)]
pub struct VerifyInput {
  trace_len: u32,
  transition_constraints: Vec<HashMap<Vec<u32>, BigUint>>,
  boundary: Vec<(u32, u32, BigUint)>
}

#[wasm_bindgen]
pub fn prove(input: JsValue) -> String {

  let input: ProveInput = serde_wasm_bindgen::from_value(input).unwrap();

  let p = BigInt::from(1) + BigInt::from(407) * BigInt::from(2).pow(119);
  let g = BigInt::from(85408008396924667383611388730472331217_u128);
  let f = Rc::new(Field::new(p, g.clone()));

  let stark = Stark::new(
    &g.clone(),
    &f,
    1,
    input.trace.len().try_into().unwrap(),
    4,
    8,
    2
  );

  let transition_constraints = input.transition_constraints
    .iter()
    .map(|v| MPolynomial::from_map(&v, &f))
    .collect();
  let boundary_constraints = input.boundary
    .iter()
    .map(|(v1, v2, v3)| (v1.clone(), v2.clone(), v3.clone().to_bigint().unwrap()))
    .collect();
  let trace = input.trace
    .iter()
    .map(|v| v.iter().map(|v2| v2.to_bigint().unwrap()).collect())
    .collect();

  stark.prove(&trace, &transition_constraints, &boundary_constraints)
}

#[wasm_bindgen]
pub fn verify(proof: String, input: JsValue) {
  let input: VerifyInput = serde_wasm_bindgen::from_value(input).unwrap();

  let p = BigInt::from(1) + BigInt::from(407) * BigInt::from(2).pow(119);
  let g = BigInt::from(85408008396924667383611388730472331217_u128);
  let f = Rc::new(Field::new(p, g.clone()));

  let stark = Stark::new(
    &g.clone(),
    &f,
    1,
    input.trace_len,
    4,
    8,
    2
  );

  let transition_constraints: Vec<MPolynomial> = input.transition_constraints
    .iter()
    .map(|v| MPolynomial::from_map(&v, &f))
    .collect();
  let boundary_constraints: Vec<(u32, u32, BigInt)> = input.boundary
    .iter()
    .map(|(v1, v2, v3)| (v1.clone(), v2.clone(), v3.clone().to_bigint().unwrap()))
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
