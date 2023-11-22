use num_bigint::BigInt;
use rstark::field::Field;
use rstark::mpolynomial::MPolynomial;
use rstark::stark::Stark;
use rstark::BigIntElement;
use std::rc::Rc;

fn main() {
    let p = BigIntElement(BigInt::from(1) + BigInt::from(407) * BigInt::from(2).pow(119));
    let g = BigIntElement(BigInt::from(85408008396924667383611388730472331217_u128));
    let f = Rc::new(Field::new(p, g.clone()));

    let register_count = 40;
    let sequence_len = 40;
    let stark = Stark::new(&g.clone(), &f, register_count, sequence_len, 32, 26, 2);

    let first_step = (0..register_count)
        .map(|v| BigIntElement(2 + BigInt::from(v)))
        .collect();

    let mut trace: Vec<Vec<BigIntElement>> = Vec::new();
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
    let mut transition_constraints: Vec<MPolynomial<BigIntElement>> = Vec::new();
    for i in 0..register_count_usize {
        let mut c = prev_state[i].clone();
        c.mul(&prev_state[i]);
        c.sub(&next_state[i]);
        transition_constraints.push(c);
    }

    let proof = stark.prove(&trace, &transition_constraints, &boundary_constraints);
    stark.verify(&proof, &transition_constraints, &boundary_constraints);
}
