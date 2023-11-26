use std::marker::PhantomData;

use crate::FieldElement;

pub struct Tree<T: FieldElement> {
    pub levels: Vec<Vec<[u8; 32]>>,
    data: PhantomData<T>,
}

impl<T: FieldElement> Tree<T> {
    pub fn hash(leaf1: &[u8; 32], leaf2: &[u8; 32]) -> [u8; 32] {
        let mut hasher = blake3::Hasher::new();
        hasher.update(leaf1);
        hasher.update(leaf2);
        *hasher.finalize().as_bytes()
    }

    pub fn root(&self) -> [u8; 32] {
        self.levels[self.levels.len() - 1][0]
    }

    pub fn leaves(&self) -> &Vec<[u8; 32]> {
        &self.levels[0]
    }

    pub fn build(leaves: &Vec<[u8; 32]>) -> Tree<T> {
        let mut levels = Vec::new();
        levels.push(leaves.clone());

        // zzzz
        let level_count = (levels[0].len() as f32).log2().ceil() as usize;

        for i in 0..level_count {
            let mut level = Vec::new();
            if levels[i].len() % 2 == 1 {
                levels[i].push(vec![0_u8; 32].try_into().unwrap());
            }
            for j in (0..levels[i].len()).step_by(2) {
                level.push(Tree::<T>::hash(&levels[i][j], &levels[i][j + 1]));
            }
            levels.push(level);
        }
        Tree {
            levels,
            data: PhantomData,
        }
    }

    pub fn commit_elements(leaves: &[T]) -> [u8; 32] {
        Self::commit(
            &leaves
                .iter()
                .map(|t| t.to_bytes_le_sized())
                .collect::<Vec<[u8; 32]>>(),
        )
    }

    pub fn commit(leaves: &Vec<[u8; 32]>) -> [u8; 32] {
        let tree = Self::build(leaves);
        tree.root()
    }

    pub fn open(&self, index: u32) -> (Vec<[u8; 32]>, [u8; 32]) {
        let mut index = index;
        if index > self.leaves().len().try_into().unwrap() {
            panic!("index is greater than leaves length");
        }
        let mut path = Vec::new();
        for i in 0..(self.levels.len() - 1) {
            let sibling_index = if index % 2 == 0 { index + 1 } else { index - 1 };
            let sibling = self.levels[i][usize::try_from(sibling_index).unwrap()];
            let node = self.levels[i][usize::try_from(index).unwrap()];
            if index % 2 == 0 {
                path.push(node);
                path.push(sibling);
            } else {
                path.push(sibling);
                path.push(node);
            }
            index >>= 1;
        }
        (path, self.root())
    }

    pub fn verify(root: &[u8; 32], _index: u32, path: &[[u8; 32]], leaf: &[u8; 32]) -> bool {
        let mut index = _index;
        let mut calculated_root = *leaf;
        for p in path.chunks(2) {
            let node_index = index % 2;
            if p[node_index as usize] != calculated_root {
                panic!("Invalid intermediate root");
            }
            calculated_root = Self::hash(&p[0], &p[1]);
            index >>= 1;
        }
        if &calculated_root != root {
            panic!("root mismatch");
        }
        true
    }
}

pub fn u128_to_bytes(v: &u128) -> [u8; 32] {
    let bytes = v.to_le_bytes();
    let mut bytes_vec = vec![0_u8; 16];
    bytes_vec.copy_from_slice(&bytes);
    bytes_vec.resize(32, 0);
    bytes_vec[0..32].try_into().unwrap()
}

#[cfg(test)]
mod tests {
    use crate::CryptoBigIntElement;

    use super::*;

    #[test]
    fn should_build_merkle_tree() {
        let mut leaves = Vec::new();
        for i in 0..100 {
            leaves.push(u128_to_bytes(&(i as u128)));
        }
        let tree = Tree::<CryptoBigIntElement>::build(&leaves);
        // expected levels.len() = log2(leaves.len()) + 1
        let expected_len = 7 + 1;
        assert_eq!(tree.levels.len(), expected_len);
        // intermediate levels should have even number of leaves
        for i in 0..(expected_len - 1) {
            assert_eq!(tree.levels[i].len() % 2, 0);
        }
        // top level should be root
        assert_eq!(tree.levels[expected_len - 1].len(), 1);
    }

    #[test]
    fn should_commit_root() {
        let mut leaves = Vec::new();

        for i in 0..100 {
            leaves.push(u128_to_bytes(&(i as u128)));
        }
        let root = Tree::<CryptoBigIntElement>::commit(&leaves);
        // choose a random value to change root
        leaves[0] = u128_to_bytes(&124812491);
        let root_changed = Tree::<CryptoBigIntElement>::commit(&leaves);

        assert_ne!(root, root_changed);
    }

    #[test]
    fn should_open_verify_tree() {
        let mut leaves = Vec::new();
        for i in 0..100 {
            leaves.push(u128_to_bytes(&(i as u128)));
        }
        let index = 5;
        let tree = Tree::<CryptoBigIntElement>::build(&leaves);
        let (path, root) = tree.open(index);
        Tree::<CryptoBigIntElement>::verify(&root, index, &path, &leaves[index as usize]);
    }

    #[test]
    #[should_panic]
    fn should_fail_to_verify() {
        let mut leaves = Vec::new();
        for i in 0..100 {
            leaves.push(u128_to_bytes(&(i as u128)));
        }
        let index = 5;
        let tree = Tree::<CryptoBigIntElement>::build(&leaves);
        let (mut path, root) = tree.open(index);
        // change some path element
        path[4] = u128_to_bytes(&124812491);
        Tree::<CryptoBigIntElement>::verify(&root, index, &path, &leaves[index as usize]);
    }

    #[test]
    #[should_panic]
    fn should_fail_to_verify_root() {
        let mut leaves = Vec::new();
        for i in 0..100 {
            leaves.push(u128_to_bytes(&(i as u128)));
        }
        let tree = Tree::<CryptoBigIntElement>::build(&leaves);
        let index = 5;
        let (path, _) = tree.open(index);
        let root = u128_to_bytes(&1921);
        Tree::<CryptoBigIntElement>::verify(&root, index, &path, &leaves[index as usize]);
    }
}
