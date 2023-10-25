use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct Message {
    pub data: Vec<u8>,
}

#[derive(Default)]
pub struct Channel {
    pub messages: Vec<Message>,
    read_index: usize,
}

impl Channel {
    pub fn new() -> Channel {
        Channel {
            messages: Vec::new(),
            read_index: 0,
        }
    }

    pub fn push(&mut self, message: &[[u8; 32]]) {
        let msg = Message {
            data: message.iter().copied().flatten().collect::<Vec<u8>>(),
        };
        self.messages.push(msg);
    }

    pub fn push_single(&mut self, message: &[u8; 32]) {
        self.push(&[*message])
    }

    pub fn pull(&mut self) -> &Message {
        // will panic if pulling past end of message vec
        let m = &self.messages[self.read_index];
        self.read_index += 1;
        m
    }

    pub fn pull_root(&mut self) -> [u8; 32] {
        let m = self.pull();
        m.data.clone().try_into().unwrap()
    }

    pub fn pull_path(&mut self) -> Vec<[u8; 32]> {
        let m = self.pull();
        let mut out: Vec<[u8; 32]> = Vec::new();
        for d in m.data.chunks(32) {
            out.push(d[0..].try_into().unwrap());
        }
        out
    }

    pub fn prover_hash(&self) -> [u8; 32] {
        let mut hasher = blake3::Hasher::new();
        for msg in &self.messages {
            for v in &msg.data {
                hasher.update(&v.to_le_bytes());
            }
        }
        *hasher.finalize().as_bytes()
    }

    pub fn verifier_hash(&self) -> [u8; 32] {
        let mut hasher = blake3::Hasher::new();
        for msg in &self.messages[0..self.read_index] {
            for v in &msg.data {
                hasher.update(&v.to_le_bytes());
            }
        }
        *hasher.finalize().as_bytes()
    }

    pub fn serialize(&self) -> String {
        serde_json::to_string(&self.messages).unwrap()
    }

    pub fn deserialize(data: &str) -> Channel {
        Channel {
            messages: serde_json::from_str(data).unwrap(),
            read_index: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::tree::u128_to_bytes;

    use super::*;

    #[test]
    fn should_generate_verifier_hash() {
        let mut c = Channel::new();
        let mut test_vec1 = Vec::new();
        let mut test_vec2 = Vec::new();
        for x in 0..20 {
            test_vec1.push(u128_to_bytes(&x));
            test_vec2.push(u128_to_bytes(&x));
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
        let mut test_vec = Vec::new();
        for x in 0..20 {
            test_vec.push(u128_to_bytes(&x));
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
        let mut test_vec = Vec::new();
        for x in 0..32 {
            test_vec.push(x as u8);
        }
        c.push_single(&test_vec.clone().try_into().unwrap());
        c.push_single(&test_vec.clone().try_into().unwrap());
        let out = &c.pull().data;
        assert_eq!(out.len(), test_vec.len());
        for x in 0..out.len() {
            assert_eq!(out[x], test_vec[x]);
        }
    }

    #[test]
    #[should_panic]
    fn should_fail_to_pull() {
        let mut c = Channel::new();
        let test_vec = [0 as u8; 32];
        c.push_single(&test_vec);
        // pull the message
        c.pull();
        // try to pull another
        c.pull();
    }
}
