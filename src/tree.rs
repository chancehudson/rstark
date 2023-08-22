use num_bigint::BigUint;

pub struct Tree {
  pub levels: Vec<Vec<BigUint>>
}

impl Tree {
  pub fn hash(leaf1: &BigUint, leaf2: &BigUint) -> BigUint {
    let mut hasher = blake3::Hasher::new();
    hasher.update(&leaf1.to_bytes_le());
    hasher.update(&leaf2.to_bytes_le());
    BigUint::from_bytes_le(hasher.finalize().as_bytes())
  }

  pub fn root(&self) -> BigUint {
    self.levels[self.levels.len() - 1][0].clone()
  }

  pub fn leaves(&self) -> &Vec<BigUint> {
    &self.levels[0]
  }

  pub fn build(leaves: &Vec<BigUint>) -> Tree {
    let mut levels: Vec<Vec<BigUint>> = Vec::new();
    levels.push(leaves.clone());

    // zzzz
    let level_count = (levels[0].len() as f32).log2().ceil() as usize;

    for i in 0..level_count {
      let mut level: Vec<BigUint> = Vec::new();
      if levels[i].len() % 2 == 1 {
        levels[i].push(BigUint::from(0 as u32));
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

  pub fn commit(leaves: &Vec<BigUint>) -> BigUint {
    let tree = Self::build(leaves);
    tree.root()
  }

  pub fn open(&self, index: u32) -> (Vec<BigUint>, BigUint) {
    let mut index = index;
    if index > self.leaves().len().try_into().unwrap() {
      panic!("index is greater than leaves length");
    }
    let mut path: Vec<BigUint> = Vec::new();
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

  pub fn verify(root: &BigUint, _index: u32, path: &Vec<BigUint>, leaf: &BigUint) -> bool {
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
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn should_build_merkle_tree() {
    let mut leaves: Vec<BigUint> = Vec::new();
    for i in 0..100 {
      leaves.push(BigUint::from(i as u32));
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
    let mut leaves: Vec<BigUint> = Vec::new();

    for i in 0..100 {
      leaves.push(BigUint::from(i as u32));
    }
    let root = Tree::commit(&leaves);
    // choose a random value to change root
    leaves[0] = BigUint::from(21940124 as u32);
    let root_changed = Tree::commit(&leaves);

    assert_ne!(root, root_changed);
  }

  #[test]
  fn should_open_verify_tree() {
    let mut leaves: Vec<BigUint> = Vec::new();
    for i in 0..100 {
      leaves.push(BigUint::from(i as u32));
    }
    let index = 5;
    let tree = Tree::build(&leaves);
    let (path, root) = tree.open(index);
    Tree::verify(&root, index, &path, &leaves[index as usize]);
  }

  #[test]
  #[should_panic]
  fn should_fail_to_verify() {
    let mut leaves: Vec<BigUint> = Vec::new();
    for i in 0..100 {
      leaves.push(BigUint::from(i as u32));
    }
    let index = 5;
    let tree = Tree::build(&leaves);
    let (mut path, root) = tree.open(index);
    // change some path element
    path[4] = BigUint::from(1290 as u32);
    Tree::verify(&root, index, &path, &leaves[index as usize]);
  }

  #[test]
  #[should_panic]
  fn should_fail_to_verify_root() {
    let mut leaves: Vec<BigUint> = Vec::new();
    for i in 0..100 {
      leaves.push(BigUint::from(i as u32));
    }
    let tree = Tree::build(&leaves);
    let index = 5;
    let (path, mut root) = tree.open(index);
    root += BigUint::from(1 as u32);
    Tree::verify(&root, index, &path, &leaves[index as usize]);
  }
}
