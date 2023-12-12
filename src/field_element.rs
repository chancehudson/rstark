use crypto_bigint::modular::runtime_mod::{DynResidue, DynResidueParams};
use crypto_bigint::{Encoding, U256, U128, U64};
use serde::{Deserialize, Serialize};
use std::hash::Hash;
use std::fmt::Debug;

/// 64 bit field
// pub const P: ParamWrapper = ParamWrapper(DynResidueParams::new(&UC::from_u64(18446744069414584321_u64)));
// pub const G: CryptoBigIntElement = CryptoBigIntElement(DynResidue::new(&UC::from_u64(2717_u64), P.0));

/// 120 bit field
pub const P: ParamWrapper = ParamWrapper(DynResidueParams::new(&UC::from_u128(1_u128 + 407_u128 * 2_u128.pow(119))));
pub const G: CryptoBigIntElement = CryptoBigIntElement(DynResidue::new(&UC::from_u128(85408008396924667383611388730472331217_u128), P.0));

pub type UC = U128;
pub const P_BITS: usize = 128;

#[cfg(target_pointer_width = "16")]
const POINTER_WIDTH: usize = 16;

#[cfg(target_pointer_width = "32")]
const POINTER_WIDTH: usize = 32;

#[cfg(target_pointer_width = "64")]
const POINTER_WIDTH: usize = 64;

pub const LIMBS: usize = P_BITS / POINTER_WIDTH;

pub trait FieldElement: Eq + PartialEq + Clone + PartialOrd + Hash + Debug {
    type ParamsType: Serialize + Debug;
    fn add(&self, v: &Self) -> Self;
    fn sub(&self, v: &Self) -> Self;
    fn mul(&self, v: &Self) -> Self;
    fn div(&self, v: &Self) -> Self;
    fn modpow(&self, e: &Self) -> Self;
    fn inv(&self) -> Self;
    fn neg(&self) -> Self;

    fn zero(p: &Self::ParamsType) -> Self;
    fn one(p: &Self::ParamsType) -> Self;
    fn two(p: &Self::ParamsType) -> Self;

    fn bits(&self) -> u64;
    fn from_u32(v: u32, p: &Self::ParamsType) -> Self;
    fn from_i32(v: i32, p: &Self::ParamsType) -> Self;
    fn to_u32(&self) -> u32;
    fn from_bytes_le(v: &[u8], p: &Self::ParamsType) -> Self;
    fn to_bytes_le(&self) -> Vec<u8>;
    fn to_bytes_le_sized(&self) -> [u8; 32];

    fn from_params(v: &Self::ParamsType) -> Self;
    fn get_params(&self) -> Self::ParamsType;
}


#[derive(Clone, PartialEq, Eq, Debug)]
pub struct CryptoBigIntElement(pub DynResidue<LIMBS>);

impl<'de> Deserialize<'de> for CryptoBigIntElement {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::de::Deserializer<'de>,
    {
        let mut bytes = Vec::deserialize(deserializer)?;
        bytes.resize(32, 0);
        Ok(CryptoBigIntElement::from_bytes_le(
            &bytes[..],
            &P
        ))
    }
}

impl PartialOrd for CryptoBigIntElement {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.0.retrieve().cmp(&other.0.retrieve()))
    }
}

impl Hash for CryptoBigIntElement {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.retrieve().hash(state)
    }
}

impl FieldElement for CryptoBigIntElement {
    type ParamsType = ParamWrapper;

    fn add(&self, v: &Self) -> Self {
        CryptoBigIntElement(self.0 + v.0)
    }

    fn sub(&self, v: &Self) -> Self {
        CryptoBigIntElement(self.0 - v.0)
    }

    fn mul(&self, v: &Self) -> Self {
        CryptoBigIntElement(self.0 * v.0)
    }

    fn div(&self, v: &Self) -> Self {
        CryptoBigIntElement(self.0 * v.0.invert().0)
    }

    fn modpow(&self, e: &Self) -> Self {
        CryptoBigIntElement(self.0.pow(&e.0.retrieve()))
    }

    fn inv(&self) -> Self {
        CryptoBigIntElement(self.0.invert().0)
    }

    fn neg(&self) -> Self {
        CryptoBigIntElement(self.0.neg())
    }

    fn zero(p: &ParamWrapper) -> Self {
        CryptoBigIntElement(DynResidue::new(&UC::ZERO, p.0))
    }

    fn one(p: &ParamWrapper) -> Self {
        CryptoBigIntElement(DynResidue::new(&UC::ONE, p.0))
    }

    fn two(p: &ParamWrapper) -> Self {
        CryptoBigIntElement(DynResidue::new(&UC::from_u32(2), p.0))
    }

    fn to_u32(&self) -> u32 {
        u32::from_le_bytes(self.to_bytes_le_sized()[..LIMBS].try_into().unwrap())
    }

    fn from_bytes_le(v: &[u8], p: &ParamWrapper) -> Self {
        // println!("{:?}", v.len());
        CryptoBigIntElement(DynResidue::new(&UC::from_le_slice(&v[0..(P_BITS/8)]), p.0))
    }

    fn to_bytes_le(&self) -> Vec<u8> {
        self.0.retrieve().to_le_bytes().to_vec()
    }

    fn to_bytes_le_sized(&self) -> [u8; 32] {
        let mut extended = self.to_bytes_le();
        extended.resize(32, 0);
        extended.try_into().unwrap()
    }

    fn from_i32(v: i32, p: &ParamWrapper) -> Self {
        if v.is_negative() {
            let s = DynResidue::new(&UC::from_u32(-v as u32), p.0);
            // let s = -v
            // calculate (s/p + 1)*p - s
            CryptoBigIntElement(
                (s * Self::from_params(p).inv().0 + Self::one(p).0) * Self::from_params(p).0 - s,
            )
        } else {
            let val = DynResidue::new(&UC::from_u32(v as u32), p.0);
            CryptoBigIntElement(val)
        }
    }

    fn from_u32(v: u32, p: &ParamWrapper) -> Self {
        CryptoBigIntElement(DynResidue::new(&UC::from_u32(v), p.0))
    }

    fn bits(&self) -> u64 {
        self.0.retrieve().bits().try_into().unwrap()
    }

    fn from_params(v: &Self::ParamsType) -> Self {
        CryptoBigIntElement(DynResidue::new(v.0.modulus(), v.0))
    }

    fn get_params(&self) -> Self::ParamsType {
        ParamWrapper(*self.0.params())
    }
}

#[derive(Debug)]
pub struct ParamWrapper(pub DynResidueParams<LIMBS>);

impl Serialize for ParamWrapper {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.0.modulus().to_string())
    }
}
