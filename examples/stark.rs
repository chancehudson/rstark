use rstark::field::Field;
use rstark::stark::Stark;
use rstark::mpolynomial::MPolynomial;
use std::rc::Rc;
use num_bigint::{BigInt, BigUint, ToBigInt};
use serde::{Serialize, Deserialize};
use std::collections::HashMap;

fn main() {
  let p = BigInt::from(1) + BigInt::from(407) * BigInt::from(2).pow(119);
  let g = BigInt::from(85408008396924667383611388730472331217_u128);
  let f = Rc::new(Field::new(p, g.clone()));

  let sequence_len = 40;
  let stark = Stark::new(
    &g.clone(),
    &f,
    2,
    sequence_len,
    32,
    32,
    2
  );

  let mut trace: Vec<Vec<BigInt>> = Vec::new();
  trace.push(vec!(BigInt::from(2), BigInt::from(3)));
  trace.push(vec!(BigInt::from(4), BigInt::from(9)));
  while trace.len() < sequence_len.try_into().unwrap() {
    let e1 = &trace[trace.len() - 1][0];
    let e2 = &trace[trace.len() - 1][1];
    trace.push(vec!(f.mul(e1, e1), f.mul(e2, e2)));
  }

  let boundary_constraints = vec!(
    (0, 0, BigInt::from(2)),
    (0, 1, BigInt::from(3)),
    (sequence_len-1, 0, trace[trace.len()-1][0].clone()),
    (sequence_len-1, 1, trace[trace.len()-1][1].clone())
  );

  let variables = MPolynomial::variables(1+2*2, &f);

  let cycle_index = &variables[0];
  let prev_state = &variables[1..3];
  let next_state = &variables[3..];
  let mut transition_constraints: Vec<MPolynomial> = Vec::new();
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

  let proof = stark.prove(&trace, &transition_constraints, &boundary_constraints);
  stark.verify(&proof, &transition_constraints, &boundary_constraints);
}
