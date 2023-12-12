use rstark::field::Field;
use rstark::mpolynomial::MPolynomial;
use rstark::stark::Stark;
use rstark::field_element::{G, CryptoBigIntElement};
use std::rc::Rc;
use std::time::Instant;

fn main() {
    let f = Rc::new(Field::new(G));

    let register_count = 40;
    let sequence_len = 40;
    let stark = Stark::new(&G, &f, register_count, sequence_len, 32, 26, 2);

    let first_step: Vec<CryptoBigIntElement> =
        (0..register_count).map(|v| f.biguint(2 + v)).collect();

    let mut trace = Vec::new();
    trace.push(first_step);
    while trace.len() < sequence_len.try_into().unwrap() {
        let last = &trace[trace.len() - 1];
        let mut next = Vec::new();
        for i in last {
            next.push(f.mul(&i, &i));
        }
        trace.push(next);
    }
    let mut boundary_constraints = Vec::new();
    for (i, v) in trace[0].iter().enumerate() {
        boundary_constraints.push((0, i as u32, v.clone()));
    }
    for (i, v) in trace[trace.len() - 1].iter().enumerate() {
        boundary_constraints.push((register_count - 1, i as u32, v.clone()));
    }

    let variables = MPolynomial::variables(1 + 2 * register_count, &f);

    let register_count_usize = usize::try_from(register_count).unwrap();

    let _cycle_index = &variables[0];
    let prev_state = &variables[1..(1 + register_count_usize)];
    let next_state = &variables[(1 + register_count_usize)..];
    let mut transition_constraints = Vec::new();
    for i in 0..register_count_usize {
        let mut c = prev_state[i].clone();
        c.mul(&prev_state[i]);
        c.sub(&next_state[i]);
        transition_constraints.push(c);
    }

    let now = Instant::now();

    let proof = stark.prove(&trace, &transition_constraints, &boundary_constraints);

    println!("Proving time: {:.2?}", now.elapsed());

    stark.verify(&proof, &transition_constraints, &boundary_constraints);
}
