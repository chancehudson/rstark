use num_bigint::BigInt;
use num_bigint::Sign;

pub struct Tree {}

impl Tree {
  pub fn hash(leaf1: &BigInt, leaf2: &BigInt) -> BigInt {
    let mut hasher = blake3::Hasher::new();
    hasher.update(&leaf1.to_bytes_le().1);
    hasher.update(&leaf2.to_bytes_le().1);
    BigInt::from_bytes_le(Sign::Plus, hasher.finalize().as_bytes())
  }

  pub fn build(leaves: &Vec<BigInt>) -> Vec<Vec<BigInt>> {

    let mut levels: Vec<Vec<BigInt>> = Vec::new();
    levels.push(leaves.clone());

    // zzzz
    let level_count = (leaves.len() as f32).log2().ceil() as usize;

    for i in 0..level_count {
      let mut level: Vec<BigInt> = Vec::new();
      if levels[i].len() % 2 == 1 {
        levels[i].push(BigInt::from(0));
      }
      for j in (0..levels[i].len()).step_by(2) {
        level.push(Tree::hash(&levels[i][j], &levels[i][j+1]));
      }
      levels.push(level);
    }
    levels
  }

  pub fn commit(leaves: &Vec<BigInt>) -> BigInt {
    let levels = Self::build(leaves);
    if levels.len() < 1 {
      panic!("invalid tree height");
    }
    levels[levels.len() - 1][0].clone()
  }

  pub fn open(_index: usize, leaves: &Vec<BigInt>) -> (Vec<Vec<BigInt>>, BigInt) {
    let tree = Self::build(leaves);
    let mut index = _index;
    if index > leaves.len() {
      panic!("index is greater than leaves length");
    }
    let mut path: Vec<Vec<BigInt>> = Vec::new();
    for i in 0..(tree.len() - 1) {
      let sibling_index;
      if index % 2 == 0 {
        sibling_index = index + 1;
      } else {
        sibling_index = index - 1;
      }
      let sibling = tree[i][sibling_index].clone();
      let node = tree[i][index].clone();
      if index % 2 == 0 {
        path.push(vec!(node, sibling));
      } else {
        path.push(vec!(sibling, node));
      }
      index >>= 1;
    }
    (path, tree[tree.len() - 1][0].clone())
  }

  pub fn verify(root: &BigInt, _index: usize, path: &Vec<Vec<BigInt>>, leaf: &BigInt) -> bool {
    let mut index = _index;
    let mut calculated_root = leaf.clone();
    for p in path {
      let node_index = index % 2;
      if p[node_index] != calculated_root {
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
    let mut leaves: Vec<BigInt> = Vec::new();
    for i in 0..100 {
      leaves.push(BigInt::from(i));
    }
    let levels = Tree::build(&leaves);
    // expected levels.len() = log2(leaves.len()) + 1
    let expected_len = 7 + 1;
    assert_eq!(levels.len(), expected_len);
    // intermediate levels should have even number of leaves
    for i in 0..(expected_len-1) {
      assert_eq!(levels[i].len() % 2, 0);
    }
    // top level should be root
    assert_eq!(levels[expected_len-1].len(), 1);
  }

  #[test]
  fn should_commit_root() {
    let mut leaves: Vec<BigInt> = Vec::new();

    for i in 0..100 {
      leaves.push(BigInt::from(i));
    }
    let root = Tree::commit(&leaves);
    // choose a random value to change root
    leaves[0] = BigInt::from(21940124);
    let root_changed = Tree::commit(&leaves);

    assert_ne!(root, root_changed);
  }

  #[test]
  fn should_open_verify_tree() {
    let mut leaves: Vec<BigInt> = Vec::new();
    for i in 0..100 {
      leaves.push(BigInt::from(i));
    }
    let index = 5;
    let (path, root) = Tree::open(index, &leaves);
    Tree::verify(&root, index, &path, &leaves[index]);
  }

  #[test]
  #[should_panic]
  fn should_fail_to_verify() {
    let mut leaves: Vec<BigInt> = Vec::new();
    for i in 0..100 {
      leaves.push(BigInt::from(i));
    }
    let index = 5;
    let (mut path, root) = Tree::open(index, &leaves);
    // change some path element
    path[4][0] = BigInt::from(1290);
    Tree::verify(&root, index, &path, &leaves[index]);
  }

  #[test]
  #[should_panic]
  fn should_fail_to_verify_root() {
    let mut leaves: Vec<BigInt> = Vec::new();
    for i in 0..100 {
      leaves.push(BigInt::from(i));
    }
    let index = 5;
    let (path, mut root) = Tree::open(index, &leaves);
    root += BigInt::from(1);
    Tree::verify(&root, index, &path, &leaves[index]);
  }
}
