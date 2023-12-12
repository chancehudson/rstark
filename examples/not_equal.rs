use rstark::field::Field;
use rstark::mpolynomial::MPolynomial;
use rstark::stark::Stark;
use rstark::field_element::{CryptoBigIntElement, ParamWrapper, UC};
use std::rc::Rc;
use crypto_bigint::modular::runtime_mod::{DynResidue, DynResidueParams};

fn main() {
    let p = ParamWrapper(DynResidueParams::new(&UC::from_u128(1_u128 + 407_u128 * 2_u128.pow(119))));
    let g = CryptoBigIntElement(DynResidue::new(&UC::from_u128(85408008396924667383611388730472331217_u128), p.0));
    let f = Rc::new(Field::new(g.clone()));
    let register_count = 3;
    let sequence_len = 2;
    let stark = Stark::new(&g.clone(), &f, register_count, sequence_len, 128, 18, 2);
    // generate a proof that two numbers (a, b) are not equal
    //
    // to do this use 3 registers
    // the first two registers store a and b
    // the next register stores 1/(a-b), let's call this z
    // the values are constrained such that
    // 0 = 1 - z*(a-b)
    //
    // if a and b are equal then a-b=0 and no value z
    // should exist such that 0 = 1 - z*0

    let mut trace = Vec::new();
    for _ in 0..2 {
        let a = f.random();
        let b = f.random();
        trace.push(vec![a.clone(), b.clone(), f.inv(&f.sub(&a, &b))]);
    }
    let mut boundary_constraints = Vec::new();
    for (i, v) in trace.iter().enumerate() {
        boundary_constraints.push((i as u32, 0, v[0].clone()));
    }

    let variables = MPolynomial::variables(1 + 2 * register_count, &f);

    let register_count_usize = usize::try_from(register_count).unwrap();

    let _cycle_index = &variables[0];
    let prev_state = &variables[1..(1 + register_count_usize)];
    let _next_state = &variables[(1 + register_count_usize)..];
    let mut transition_constraints = Vec::new();
    {
        let mut one = MPolynomial::new(&f);
        one.term(&f.one(), &vec![0]);
        one.sub(
            &prev_state[2]
                .clone()
                .mul(&prev_state[0].clone().sub(&prev_state[1])),
        );
        transition_constraints.push(one);
    }

    println!("building proof...");
    let proof = stark.prove(&trace, &transition_constraints, &boundary_constraints);
    println!("verifying proof...");
    stark.verify(&proof, &transition_constraints, &boundary_constraints);
    println!("proof valid!");
}
