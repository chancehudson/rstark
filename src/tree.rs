

pub struct Tree {
  pub levels: Vec<Vec<[u8; 32]>>
}

impl Tree {
  pub fn hash(leaf1: &[u8; 32], leaf2: &[u8; 32]) -> [u8; 32] {
    let mut hasher = blake3::Hasher::new();
    hasher.update(leaf1);
    hasher.update(leaf2);
    hasher.finalize().as_bytes().clone()
  }

  pub fn root(&self) -> [u8; 32] {
    self.levels[self.levels.len() - 1][0].clone()
  }

  pub fn leaves(&self) -> &Vec<[u8; 32]> {
    &self.levels[0]
  }

  pub fn build_elements(leaves: &Vec<u128>) -> Tree {
    Self::build(&Self::u128s_to_bytes(leaves))
  }

  pub fn build(leaves: &Vec<[u8; 32]>) -> Tree {
    let mut levels: Vec<Vec<[u8; 32]>> = Vec::new();
    levels.push(leaves.clone());

    // zzzz
    let level_count = (levels[0].len() as f32).log2().ceil() as usize;

    for i in 0..level_count {
      let mut level: Vec<[u8; 32]> = Vec::new();
      if levels[i].len() % 2 == 1 {
        levels[i].push(vec!(0_u8; 32).try_into().unwrap());
      }
      for j in (0..levels[i].len()).step_by(2) {
        level.push(Tree::hash(&levels[i][j], &levels[i][j+1]));
      }
      levels.push(level);
    }
    Tree {
      levels
    }
  }

  pub fn commit_elements(leaves: &Vec<u128>) -> [u8; 32] {
    Self::commit(&Self::u128s_to_bytes(leaves))
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
    let mut path: Vec<[u8; 32]> = Vec::new();
    for i in 0..(self.levels.len() - 1) {
      let sibling_index;
      if index % 2 == 0 {
        sibling_index = index + 1;
      } else {
        sibling_index = index - 1;
      }
      let sibling = self.levels[i][usize::try_from(sibling_index).unwrap()].clone();
      let node = self.levels[i][usize::try_from(index).unwrap()].clone();
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

  pub fn verify(root: &[u8; 32], _index: u32, path: &Vec<[u8; 32]>, leaf: &[u8; 32]) -> bool {
    let mut index = _index;
    let mut calculated_root = leaf.clone();
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
    return true;
  }

  pub fn u128s_to_bytes(v: &Vec<u128>) -> Vec<[u8; 32]> {
    let mut out: Vec<[u8; 32]> = Vec::new();
    for i in v {
      out.push(Self::u128_to_bytes(i));
    }
    out
  }

  pub fn u128_to_bytes(v: &u128) -> [u8; 32] {
    let bytes = v.to_le_bytes();
    let mut bytes_vec = vec!(0_u8; 16);
    bytes_vec.copy_from_slice(&bytes);
    bytes_vec.resize(32, 0);
    bytes_vec[0..32].try_into().unwrap()
  }

  pub fn u32_to_bytes(v: &u32) -> [u8; 32] {
    let bytes = v.to_le_bytes();
    let mut bytes_vec = vec!(0_u8; 8);
    bytes_vec.copy_from_slice(&bytes);
    bytes_vec.resize(32, 0);
    bytes_vec[0..32].try_into().unwrap()
  }

}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn should_build_merkle_tree() {
    let mut leaves: Vec<[u8; 32]> = Vec::new();
    for i in 0..100 {
      leaves.push(Tree::u128_to_bytes(&(i as u128)));
    }
    let tree = Tree::build(&leaves);
    // expected levels.len() = log2(leaves.len()) + 1
    let expected_len = 7 + 1;
    assert_eq!(tree.levels.len(), expected_len);
    // intermediate levels should have even number of leaves
    for i in 0..(expected_len-1) {
      assert_eq!(tree.levels[i].len() % 2, 0);
    }
    // top level should be root
    assert_eq!(tree.levels[expected_len-1].len(), 1);
  }

  #[test]
  fn should_commit_root() {
    let mut leaves: Vec<[u8; 32]> = Vec::new();

    for i in 0..100 {
      leaves.push(Tree::u128_to_bytes(&(i as u128)));
    }
    let root = Tree::commit(&leaves);
    // choose a random value to change root
    leaves[0] = Tree::u128_to_bytes(&124812491);
    let root_changed = Tree::commit(&leaves);

    assert_ne!(root, root_changed);
  }

  #[test]
  fn should_open_verify_tree() {
    let mut leaves: Vec<[u8; 32]> = Vec::new();
    for i in 0..100 {
      leaves.push(Tree::u128_to_bytes(&(i as u128)));
    }
    let index = 5;
    let tree = Tree::build(&leaves);
    let (path, root) = tree.open(index);
    Tree::verify(&root, index, &path, &leaves[index as usize]);
  }

  #[test]
  #[should_panic]
  fn should_fail_to_verify() {
    let mut leaves: Vec<[u8; 32]> = Vec::new();
    for i in 0..100 {
      leaves.push(Tree::u128_to_bytes(&(i as u128)));
    }
    let index = 5;
    let tree = Tree::build(&leaves);
    let (mut path, root) = tree.open(index);
    // change some path element
    path[4] = Tree::u128_to_bytes(&124812491);
    Tree::verify(&root, index, &path, &leaves[index as usize]);
  }

  #[test]
  #[should_panic]
  fn should_fail_to_verify_root() {
    let mut leaves: Vec<[u8; 32]> = Vec::new();
    for i in 0..100 {
      leaves.push(Tree::u128_to_bytes(&(i as u128)));
    }
    let tree = Tree::build(&leaves);
    let index = 5;
    let (path, mut root) = tree.open(index);
    root = Tree::u128_to_bytes(&1921);
    Tree::verify(&root, index, &path, &leaves[index as usize]);
  }
}
