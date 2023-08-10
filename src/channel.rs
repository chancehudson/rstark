use num_bigint::BigUint;

pub struct Message {
  pub data: Vec<BigUint>
}

pub struct Channel {
  pub messages: Vec<Message>,
  read_index: usize
}

impl Channel {
  pub fn new() -> Channel {
    Channel {
      messages: Vec::new(),
      read_index: 0
    }
  }

  pub fn push(&mut self, message: &Vec<BigUint>) {
    let msg = Message {
      data: message.to_vec()
    };
    self.messages.push(msg);
  }

  pub fn pull(&mut self) -> &Message {
    // will panic if pulling past end of message vec
    let m = &self.messages[self.read_index];
    self.read_index += 1;
    m
  }

  pub fn prover_hash(&self) -> BigUint {
    let mut hasher = blake3::Hasher::new();
    for msg in &self.messages {
      for v in &msg.data {
        hasher.update(&v.to_bytes_le());
      }
    }
    BigUint::from_bytes_le(hasher.finalize().as_bytes())
  }

  pub fn verifier_hash(&self) -> BigUint {
    let mut hasher = blake3::Hasher::new();
    for msg in &self.messages[0..self.read_index] {
      for v in &msg.data {
        hasher.update(&v.to_bytes_le());
      }
    }
    BigUint::from_bytes_le(hasher.finalize().as_bytes())
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn should_generate_verifier_hash() {
    let mut c = Channel::new();
    let mut test_vec1: Vec<BigUint> = Vec::new();
    let mut test_vec2: Vec<BigUint> = Vec::new();
    for x in 0..20 {
      test_vec1.push(BigUint::new(vec!(x)));
      test_vec2.push(BigUint::new(vec!(x)));
    }
    c.push(&test_vec1);
    c.push(&test_vec2);
    let first_hash = c.verifier_hash();
    assert_ne!(c.prover_hash(), first_hash);
    c.pull();
    let next_hash = c.verifier_hash();
    assert_ne!(c.prover_hash(), next_hash);
    assert_ne!(first_hash, next_hash);
    c.pull();
    let final_hash = c.verifier_hash();
    assert_ne!(final_hash, next_hash);
    assert_eq!(c.prover_hash(), final_hash);
  }

  #[test]
  fn should_generate_prover_hash() {
    let mut c = Channel::new();
    let mut test_vec: Vec<BigUint> = Vec::new();
    for x in 0..20 {
      test_vec.push(BigUint::new(vec!(x)));
    }
    c.push(&test_vec);
    let start = c.prover_hash();
    c.pull();
    let next = c.prover_hash();
    assert_eq!(start, next);
  }

  #[test]
  fn should_push_pull_message() {
    let mut c = Channel::new();
    let mut test_vec: Vec<BigUint> = Vec::new();
    for x in 0..20 {
      test_vec.push(BigUint::new(vec!(x)));
    }
    c.push(&test_vec);
    let out = c.pull();
    assert_eq!(out.data.len(), test_vec.len());
    for x in 0..out.data.len() {
      assert_eq!(out.data[x], test_vec[x]);
    }
  }

  #[test]
  #[should_panic]
  fn should_fail_to_pull() {
    let mut c = Channel::new();
    let test_vec: Vec<BigUint> = vec!(BigUint::new(vec!()));
    c.push(&test_vec);
    // pull the message
    c.pull();
    // try to pull another
    c.pull();
  }
}
