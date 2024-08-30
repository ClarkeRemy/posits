//! This implementation of Posits follows the "Standard for Posit™ Arithmetic (2022)*"
//!
//! A primary goal is to allow the ability to run posits at compile time (Rust const context).
//! A secondary goal is to integrate this with standard library traits and portable SIMD.
//!
//! A best effort is done to follow the standard, and as a consequence will heavily quote it for clarity in the documentation.
//! A copy of the standary is included in the repository "Standard_for_Posit(TM)_Arithmetic_(2022).pdf" 
//! taken from the [site](https://posithub.org/docs/posit_standard-2.pdf).
#![allow(unused)]
#![no_implicit_prelude]
pub(crate) mod prelude;
use prelude::*;

mod posit_ops;


mod posit_bits { pub trait Sealed {} }

macro_rules! posit_internals 
{($(| $I_TYPE:ty ; $U_TYPE:ty , $U_WIDE_TYPE:ty)*) => {$(  
  impl posit_bits::Sealed for $I_TYPE {}
  impl PositBits for $I_TYPE {}
  impl PositRepr for Posit<$I_TYPE> {type IRepr = $I_TYPE; type URepr = $U_TYPE; type UWide = $U_WIDE_TYPE;}
)*};}
posit_internals!
{
| i8    ; u8    , u16 
| i16   ; u16   , u32 
| i32   ; u32   , u64 
| i64   ; u64   , u128 
| isize ; usize , [usize;2]
| i128  ; u128  , [u128 ;2]
}

/// This trait is [sealed]. It is only implemented by primitive types that are signed integers. 
/// The choice signed integer types only is because Posits are also signed, and have a total order 
/// isomorphic to the signed integer types with the same bit pattern.
///
#[doc = "[sealed]: <https://rust-lang.github.io/api-guidelines/future-proofing.html\
#sealed-traits-protect-against-downstream-implementations-c-sealed>"]
pub trait PositBits : posit_bits::Sealed {}

mod posit_repr 
{ pub trait Sealed {}
  impl<T : crate::PositBits> Sealed for crate::Posit<T> {}
}

/// Describes the types used for implementing the [`Posit`] type
pub trait PositRepr : posit_repr::Sealed
{ /** signed        repr */ type IRepr : PositBits + Copy + Eq + Debug
; /** unsigned      repr */ type URepr : Copy + Eq + Debug
; /** unsigned wide repr */ type UWide : Copy + Eq + Debug
;}

// Posits are first implemented as genericaly as possible using isize. 
// The other types are then implemented with a templated version of the same code.
#[derive(PartialEq, Eq, Clone, Copy)]
pub struct Posit<T : PositBits>(pub(crate) T);

/// used in return type of `posit.quotient_remainder(posit)`
///
/// this should only be used as a return value.
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub struct QuotRem<T : PositBits> where Posit<T> : Debug
{ pub(crate) quot : Posit<T>, pub(crate) rem : Posit<T> }

impl<T : PositBits + Copy> QuotRem<T> where Posit<T> : Debug
{ fn quotient(&self) -> Posit<T> {self.quot}
  fn remainder(&self) -> Posit<T> {self.rem}
}


/// - standard : "exponent _______ The power-of-two scaling determined by the exponent bits, in the set {0, 1, 2, 3}."
/// - standard : "exponent bits __ A two-bit unsigned integer bit field that determines the exponent."
#[repr(transparent)] // repr i8 makes calculations with the regime which is also i8 more straightforward
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub struct PositExponent(i8);

/// Posit Value exploded to its constitent parts, is never NaR or 0
///
/// most of the fields are self evident, with the exception of the fraction 
/// The representation of the fraction, the highest bit is reserved with a 0 bit; the rest of the fraction starts at the n-1 bit
/// 
/// Effectively it is a fixed point number with the first bit set as the integral part, the fraction follows.
/// 
/// example for an Posit<i8> :
/// - <Posit<i8> as PositRepr>::URepr == u8
/// - imagine a binary fixed point literal where `0b_0_0000000 => 0b_0.0...`
/// - `1/2 => 0b_0_1000000 => 0b_0.10...`
/// - `1   => 0b_1_0000000 => 0b_1.0...`
/// - given a __POSITIVE__ Posit with : 
///   - regime   = r
///   - exponent = e
///   - fraction = f
///   - value is `2.pow(4*r+e) * (1+f)`
///   - looking at just when `f = 1/2, 1+f = 0b_1.0 + 0b_0.1 = 0b_1.1`
/// 
/// This representation makes multiplication and addition relatively straightforward to implement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PositUnpacked<T : PositBits> where Posit<T> : PositRepr
{ /// true when sign bit is on, `posit < 0`
  ///
  /// standard : "sign bit__ The MSB of a posit or quire format."
  sign_bit : bool
, /// standard : "regime __ The power-of-16 scaling determined by the regime bits. It is a signed integer."
  regime : i8 // the largest number of bits will be 128 bits for i128, but the sign bit is always there, so only a count up to 127 is required, perfect for an i8
, /// standard : "exponent __ The power-of-two scaling determined by the exponent bits, in the set {0, 1, 2, 3}."
  ///
  /// standard : "exponent bits __ A two-bit unsigned integer bit field that determines the exponent."
  exponent : PositExponent
, /// 1 + fraction, the highest bit is always 1, the fraction is always the remaining bits
  ///
  /// standard : "fraction __ The value represented by the fraction bits; 0 ≤ fraction < 1."
  significand : <Posit<T> as PositRepr>::URepr
}

impl PositUnpacked<i64>
{ /// encodes a posit value that is non-zero and non-NaR
  const fn encode(self) -> Posit<i64> 
  { type P = Posit<i64>
  ; let Self { sign_bit, regime, exponent, significand } = self
  ; core::debug_assert!(regime <   ( P::PRECISION-1) as i8 )
  ; core::debug_assert!(regime >= -((P::PRECISION-1) as i8) )
  ; core::debug_assert!(exponent.0 < P::REGIME_BUMP as i8)
  ; let exp_fract = ((exponent.0 as <P as PositRepr>::URepr) << P::PRECISION-2 | (significand << 1) >> 2 )
  ; let unsigned  = if regime < 0 
         { let shift = ((regime.abs()+1) as u32)
         ; if shift > P::PRECISION 
                {1} 
           else {(1<<P::PRECISION-1 | exp_fract >> 1) >> shift}
         } // negative exponent has leading zeros
    else { let shift = (regime as u32+3)
         ; (!0<<P::PRECISION-1-regime as u32)>>1 |  exp_fract >> /* max */ if shift > P::PRECISION {0} else {shift}
         } /* sign_bit + terminator + 1 because positive starts with one bit */

  ; let mut u_posit = Posit(unsigned as <P as PositRepr>::IRepr)
  ; if sign_bit { u_posit.negate() } else { u_posit }
  }
}

macro_rules! posit_unpacked_impl
{($($TYPE:ty)*) => {$(
  impl PositUnpacked<$TYPE>
  { const fn encode(self) -> Posit<$TYPE> 
      { type P = Posit<$TYPE>
      ; let Self { sign_bit, regime, exponent, significand } = self
    ; core::debug_assert!(regime <   ( P::PRECISION-1) as i8 )
    ; core::debug_assert!(regime >= -((P::PRECISION-1) as i8) )
    ; core::debug_assert!(exponent.0 < P::REGIME_BUMP as i8)
    ; let exp_fract = ((exponent.0 as <P as PositRepr>::URepr) << P::PRECISION-2 | (significand << 1) >> 2 )
    ; let unsigned  = if regime < 0 
           { let shift = ((regime.abs()+1) as u32)
           ; if shift >= P::PRECISION 
                  {1} 
             else {(1<<P::PRECISION-1 | exp_fract >> 1) >> shift}
           } // negative exponent has leading zeros
      else { let shift = (regime as u32+3)
           ; (!0<<P::PRECISION-1-regime as u32)>>1 |  exp_fract >> /* max */ if shift > P::PRECISION-1 {0} else {shift}
           } /* sign_bit + terminator + 1 because positive starts with one bit */
   
    ; let mut u_posit = Posit(unsigned as <P as PositRepr>::IRepr)
    ; if sign_bit { u_posit.negate() } else { u_posit }
    }
  }
)*};}
posit_unpacked_impl!{i8 i16 i32 i128 isize}


#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PositDecoded<T : PositBits> where Posit<T> : PositRepr 
{ NaR
, Zero
, NonZero(PositUnpacked<T>)
}

impl<T : PositBits> Posit<T>
{ // the power of 2 that a regime represents
  const REGIME_BUMP : u32 = 0b_1_0000_u8.trailing_zeros() // 1 regime = 2^4 = 0b_1_0000
; pub const fn from_bits(inner : T)->Self {Self(inner)}
  pub const fn to_bits(self)->T where T : Copy {self.0}

}

impl Posit<i64>
{ /// standard : "precision __ The total storage size for expressing any number format, in bits. For a posit, precision is 𝑛 bits"
  pub const PRECISION : u32  = <Self as PositRepr>::IRepr::BITS
; /// standard : "pIntMax __ The largest consecutive integer-valued posit value. It is a function of 𝑛."
  pub const P_INT_MAX : Self = Self( (2 as <Self as PositRepr>::IRepr).pow((4 * (Self::PRECISION - 3)) / 5 )  )
; /// Integer where all Posits after this one are integers.
  pub const FRACTION_LIMIT : Self = 
  { // based on the observations for the limit :
    // - `Lemma : Limit is Positive => power_of_2 = (4*regime + exponent)` 
    // - `forall posit_integers => power_of_2 >= the fraction_bit_len`.
    // - `forall exponent => (0..=3).conains(exponent)`
    //
    // the subtletly is to think iteratively :
    // - every time the regime becomes 1 bit longer, the fraction must be 1 bit shorter.
    // - every regime is worth 4 `powers of 2`, plus the removal of 1 bit. So extention of the regime negates 5 fraction bits
    // - iterative subtraction can be done with integer division `/ 5bits`
    // - the remainder `% 5bits` could be 4,
    //   but the exponent can only constitute up to 3 `powers of 2`.
    //   In the edge case of 4, the regime must be bumped, leaving the exponent 0

    let max_fract_bits = Self::PRECISION - 5 // 5 = signbit + reg_min_bits + exponent_bits


  ; let rem      =  max_fract_bits % (Self::REGIME_BUMP+1)
  ; let regime   = (max_fract_bits / (Self::REGIME_BUMP+1))  +  rem / Self::REGIME_BUMP
  ; let exponent = rem % Self::REGIME_BUMP
  ; Self
    ( (  ((0b_10_00 + exponent as <Self as PositRepr>::IRepr) << Self::PRECISION-4) 
      >> regime+1  // when positive, regime = k - 1, so invert to get k = regime+1
      ) 
    ^ (1 << Self::PRECISION-1) // zero the sign bit making
    )
  }


; /** 0b_01...1 , MAX      */ pub const MAX_POS : Self = Self( !(1 << Self::PRECISION - 1) )
; /** 0b_010..0 , +1       */ pub const ONE     : Self = Self( 0b_01 << Self::PRECISION - 2 ) 
; /** 0b_0...01 , next(0)  */ pub const MIN_POS : Self = Self( 1 )
; /** 0b_0....0 , 0        */ pub const ZERO    : Self = Self( 0 )
; /** 0b_1....1 , prior(0) */ pub const MAX_NEG : Self = Self( -1)
; /** 0b_110..0 , -1       */ pub const NEG_ONE : Self = Self( 0b_11 << Self::PRECISION - 2 )
; /** 0b_10..01 , MIN      */ pub const MIN_NEG : Self = Self( (1 << Self::PRECISION - 1) | 1)
; /// 0b_10...0 , NAR      ,  "Not a Real" : The Error value of a posit, 
                              pub const NAR     : Self = Self( 1 << Self::PRECISION - 1 )
;

  /** 0b_00111..0  , 1/2 */ pub const HALF         : Self = Self( 0b_00111 << Self::PRECISION-5 )
; /** 0b_01000_1..0, 3/2 */ pub const THREE_HALVES : Self = Self( 0b_01000_1 << Self::PRECISION-6 )
;
}

macro_rules! posit_assoc_consts 
{($($TYPE:ty)*) => {$(
  impl Posit<$TYPE>
  { /// standard : "precision __ The total storage size for expressing any number format, in bits. For a posit, precision is 𝑛 bits"
    pub const PRECISION : u32  = <Self as PositRepr>::IRepr::BITS
  ; /// standard : "pIntMax __ The largest consecutive integer-valued posit value. It is a function of 𝑛."
    pub const P_INT_MAX : Self = Self( (2 as <Self as PositRepr>::IRepr).pow((4 * (Self::PRECISION - 3)) / 5 )  )
  ; /// Integer where all Posits after this one are integers.
    pub const FRACTION_LIMIT : Self = 
    { // based on the observations for the limit :
      // - `Lemma : Limit is Positive => power_of_2 = (4*regime + exponent)` 
      // - `forall posit_integers => power_of_2 >= the fraction_bit_len`.
      // - `forall exponent => (0..=3).conains(exponent)`
      //
      // the subtletly is to think iteratively :
      // - every time the regime becomes 1 bit longer, the fraction must be 1 bit shorter.
      // - every regime is worth 4 `powers of 2`, plus the removal of 1 bit. So extention of the regime negates 5 fraction bits
      // - iterative subtraction can be done with integer division `/ 5bits`
      // - the remainder `% 5bits` could be 4,
      //   but the exponent can only constitute up to 3 `powers of 2`.
      //   In the edge case of 4, the regime must be bumped, leaving the exponent 0
  
      let max_fract_bits = Self::PRECISION - 5 // 5 = signbit + reg_min_bits + exponent_bits
  
  
    ; let rem      =  max_fract_bits % (Self::REGIME_BUMP+1)
    ; let regime   = (max_fract_bits / (Self::REGIME_BUMP+1))  +  rem / Self::REGIME_BUMP
    ; let exponent = rem % Self::REGIME_BUMP
    ; Self
      ( (  ((0b_10_00 + exponent as <Self as PositRepr>::IRepr) << Self::PRECISION-4) 
        >> regime+1  // when positive, regime = k - 1, so invert to get k = regime+1
        ) 
      ^ (1 << Self::PRECISION-1) // zero the sign bit making
      )
    }
  
  ; /** 0b_01...1 , MAX      */ pub const MAX_POS : Self = Self( !(1 << Self::PRECISION - 1) )
  ; /** 0b_010..0 , +1       */ pub const ONE     : Self = Self( 0b_01 << Self::PRECISION - 2 ) 
  ; /** 0b_0...01 , next(0)  */ pub const MIN_POS : Self = Self( 1 )
  ; /** 0b_0....0 , 0        */ pub const ZERO    : Self = Self( 0 )
  ; /** 0b_1....1 , prior(0) */ pub const MAX_NEG : Self = Self( -1)
  ; /** 0b_110..0 , -1       */ pub const NEG_ONE : Self = Self( 0b_11 << Self::PRECISION - 2 )
  ; /** 0b_10..01 , MIN      */ pub const MIN_NEG : Self = Self( (1 << Self::PRECISION - 1) | 1)
  ; /// 0b_100..0 , NAR      ,  "Not a Real" : The Error value of a posit, 
                                pub const NAR     : Self = Self( 1 << Self::PRECISION - 1 )
  ;
    /** 0b_00111..0  , 1/2 */ pub const HALF         : Self = Self( 0b_00111   << Self::PRECISION-5 )
  ; /** 0b_01000_1..0, 3/2 */ pub const THREE_HALVES : Self = Self( 0b_01000_1 << Self::PRECISION-6 )
  ;
  }

)*};}
posit_assoc_consts!{i8 i16 i32 isize i128}


impl Posit<i64>
{ ////
  // 5.2 Basic functions of one posit value argument 
  ////

  /// +x => -x, -x => +x, else unchanged
  ///
  /// standard : "negate(posit) __ returns −posit.(7)"
  ///
  /// footnote : "(7) This is the 2’s complement of the posit representation. 2’s complement affects neither 0 nor NaR, since they are unsigned."
  pub const fn negate(self : Self) -> Self 
  {  Self( self.0.wrapping_neg() ) }

  /// -x => +x, else unchanged
  ///
  /// standard : "abs(posit) __ returns negate(posit) if posit < 0, else posit"
  pub const fn abs(self : Self) -> Self { Self( self.0.wrapping_abs() ) }

  /// +x => +1 , -x => -1, else unchanged
  ///
  /// standard : "sign(posit) __ returns a posit value: 1 if posit > 0, −1 if posit < 0, or 0 if posit = 0."

  pub const fn sign(self : Self) -> Self 
  { let nar_mask = !( ((self.0 == Self::NAR.0) as <Self as PositRepr>::IRepr) << Self::PRECISION - 2 ) // NAR => 0b_101..1 , _ => 0b_1..1
  ; Self( nar_mask & (self.0.signum() << Self::PRECISION - 2)) 
  }

  /// standard : "nearestInt(posit) __ returns the integer-valued posit value nearest to posit, and returns the nearest even integer-valued posit value if two integers are equally near."
  pub const fn nearest_int(self : Self) -> Self
  { if self.0 <= Self::FRACTION_LIMIT.negate().0 || self.0 >= Self::FRACTION_LIMIT.0 { return self } // note this includes NaR
    if self.0 >  Self::THREE_HALVES.negate().0   && self.0 <  Self::HALF.negate().0  { return Self::NEG_ONE }
    if self.0 >= Self::HALF.negate().0           && self.0 <= Self::HALF.0           { return Self::ZERO }
    if self.0 >  Self::HALF.0                    && self.0 <  Self::THREE_HALVES.0   { return Self::ONE }
    
    match self.sign_regime_exponent_significand() 
    { PositDecoded::NonZero(PositUnpacked { sign_bit, regime, exponent, significand })
      =>
      { let power = regime as u32 * Self::REGIME_BUMP + exponent.0 as u32 // this could be 0
      ; core::debug_assert!{ power <= Self::PRECISION-5 }

      ; let mask = !0 >> power+2
      ; let half_bit          = 0 != significand  &  mask+1
      ; let one_bit           = 0 != significand  &  mask+1 << 1
      ; let trailing_fraction = 0 != significand  &  mask

      ; let abs = self.abs()

      ; let integer = if half_bit & (!one_bit&trailing_fraction | one_bit)
        { abs.ceil() } else { abs.floor() }

      ; if sign_bit 
        { integer.negate() } else { integer }
      }
      PositDecoded::NaR | PositDecoded::Zero => { self } // this should in practice be unreachable
    }

  }

  /// int => int , +int.fract => +int + 1, -int.fract => -int, NAR => NAR
  ///
  /// standard : "ceil(posit) __ returns the smallest integer-valued posit value greater than or equal to posit."
  pub const fn ceil(self : Self) -> Self
  { if self.0 >  Self::NEG_ONE.0                 && self.0 <= Self::ZERO.0           { return Self::ZERO }
    if self.0 >  Self::ZERO.0                    && self.0 <= Self::ONE.0            { return Self::ONE  }
    if self.0 <= Self::FRACTION_LIMIT.negate().0 || self.0 >= Self::FRACTION_LIMIT.0 { return self } // note this includes NaR
    
    match self.sign_regime_exponent_significand()
    { PositDecoded::NonZero(PositUnpacked {sign_bit, regime, exponent, significand}) // at this point, posit.abs() > 1, so regime is  positive
      => 
      { let mut power = regime as u32 * Self::REGIME_BUMP + exponent.0 as u32
      ; let mask  = !(0 as <Self as PositRepr>::URepr) >> power+1
      ; let trailing_fraction = mask & significand
      
      ; let significand = 
        if sign_bit { significand & !mask }
        else if let Some(sig) = significand.checked_add(mask) { sig &  !mask }
        else { power+=1 ; 1 << Self::PRECISION-1 }
      
      ; let regime = (power / Self::REGIME_BUMP) as i8
      ; let exponent = PositExponent((power%Self::REGIME_BUMP) as i8)
      
      ; PositUnpacked::<<Self as PositRepr>::IRepr> {sign_bit, regime, exponent, significand}.encode()
      }
      PositDecoded::NaR | PositDecoded::Zero => { self } // this should in practice be unreachable
    }
  }

  /// int => int , +int.fract => int, -int.fract => -int - 1, NAR => NAR
  ///
  /// standard : "floor(posit) __ returns the largest integer-valued posit value less than or equal to posit."
  pub const fn floor(self : Self) -> Self
  { self.negate().ceil().negate() }

  /// The inverse of `Posit::prior`
  ///
  /// standard : "next(posit) __ returns the posit value of the lexicographic successor of posit’s representation. (8)"
  ///
  /// footnote : "(8) wrapping around, if necessary"
  ///
  /// note! lexicographical order is based on the total order; with wrapping these special cases should be in mind :
  ///
  /// MAX_POS => NAR, NAR => MIN_NEG, MAX_NEG => 0, 0 => MIN_POS
  pub const fn next(self : Self) -> Self{ Self( self.0.wrapping_add(1) ) }

  /// The inverse of `Posit::next`
  /// 
  /// standard : "prior(posit) __ returns the posit value of the lexicographic predecessor of posit’s representation. (9)"
  /// - footnote : "(9) wrapping around, if necessary"
  ///
  /// note! lexicographical order is based on the total order; with wrapping these special cases should be in mind : 
  ///
  /// MIN_NEG => NAR, NAR => MAX_POS, MIN_POS => 0, 0 => MAX_NEG
  pub const fn prior(self : Self) -> Self { Self( self.0.wrapping_sub(1) ) }


  ////
  // 5.3 Comparison functions of two posit value arguments
  ////

  pub const fn compare_equal         (self : Self, rhs : Self) -> bool {/*! standard : "compareEqual(posit1,posit2)        __ returns True if posit1 = posit2, else False." */ self.0 == rhs.0 }
  pub const fn compare_not_equal     (self : Self, rhs : Self) -> bool {/*! standard : "compareNotEqual(posit1,posit2)     __ returns True if posit1 ≠ posit2, else False." */ self.0 != rhs.0 }
  pub const fn compare_greater       (self : Self, rhs : Self) -> bool {/*! standard : "compareGreater(posit1,posit2)      __ returns True if posit1 > posit2, else False." */ self.0 >  rhs.0 }
  pub const fn compare_greater_equal (self : Self, rhs : Self) -> bool {/*! standard : "compareGreaterEqual(posit1,posit2) __ returns True if posit1 ≥ posit2, else False." */ self.0 >= rhs.0 }
  pub const fn compare_less          (self : Self, rhs : Self) -> bool {/*! standard : "compareLess(posit1,posit2)         __ returns True if posit1 < posit2, else False." */ self.0 <  rhs.0 }
  pub const fn compare_less_equal    (self : Self, rhs : Self) -> bool {/*! standard : "compareLessEqual(posit1,posit2)    __ returns True if posit1 ≤ posit2, else False." */ self.0 <= rhs.0 }

  ////
  // 5.4 Aritmetic function of two posit value arguments
  ////
  
  /// standard : "addition(posit1,posit2) __ returns posit1+posit2, rounded."
  pub const fn addition(self : Self, rhs : Self) -> Self 
  { if self.0 == Self::NAR.0
    || rhs.0  == Self::NAR.0    { return Self::NAR }
  ; if self.0 == Self::ZERO.0   { return rhs}
  ; if rhs.0  == Self::ZERO.0   { return self}
  ; if self.0 == rhs.negate().0 { return Self::ZERO }

  let (big, small) = if self.0.abs() > rhs.0.abs() 
       {(self, rhs)} 
  else {(rhs, self)}    

  ; let PositDecoded::NonZero(PositUnpacked 
  { sign_bit   : sign_b
  , regime     : reg_b
  , exponent   : exp_b
  , significand : significand_b}) 
  = big.sign_regime_exponent_significand() else { core::unreachable!() }
  
  ; let PositDecoded::NonZero(PositUnpacked { 
      sign_bit   : sign_s
    , regime     : reg_s
    , exponent   : exp_s
    , significand : significand_s
    })
  = small.sign_regime_exponent_significand() else {  core::unreachable!() }

 
  // big small
  //  b+ s  =>   b+s
  // -b+-s  => -(b+s)
  //  b+-s  =>   b-s
  // -b+ s  => -(b-s)

  ; let power_b   = (reg_b as i16)*(Self::REGIME_BUMP as i16) + (exp_b.0 as i16)
  ; let power_s   = (reg_s as i16)*(Self::REGIME_BUMP as i16) + (exp_s.0 as i16)

  ; let power_diff = power_b - power_s
  ; if power_diff >= Self::PRECISION as i16 { return big }

  ; let significand_b_ = significand_b >> 1
  ; let significand_s_ = significand_s >> 1 >> power_diff

  ; let mut significand = if sign_b^sign_s
         { significand_b_ - significand_s_ }
    else { significand_b_ + significand_s_ }

  ; let leading = significand.leading_zeros()
  ; let power = power_b + 1 - leading as i16
  ; significand <<= leading

  ; let mut regime   =               (power / Self::REGIME_BUMP as i16) as i8
  ; let mut exponent = PositExponent((power % Self::REGIME_BUMP as i16) as i8)
  ; if exponent.0 < 0 { regime -=1 ; exponent.0 += Self::REGIME_BUMP as i8 }

  ; let sign_bit = sign_b

  ; PositUnpacked::<<Self as PositRepr>::IRepr> { sign_bit, regime, exponent, significand}.encode() 
  }
  /// standard : "subtraction(posit1,posit2) __ returns posit1−posit2, rounded."
  pub const fn subtraction(self : Self, rhs : Self) -> Self { self.addition(rhs.negate()) }
  /// standard : "multiplication(posit1,posit2) returns posit1×posit2, rounded."
  pub const fn multiplication(self : Self, rhs : Self) -> Self 
  { if rhs.0  == Self::NAR.0 { return Self::NAR }
  ; if rhs.0  == Self::ONE.0 { return self}
  ; if self.0 == Self::ONE.0 { return rhs}
  ; let PositDecoded::NonZero(PositUnpacked {sign_bit, regime, exponent, significand}) 
  = self.sign_regime_exponent_significand() else { return self }

  ; let PositDecoded::NonZero(PositUnpacked 
    { sign_bit   : sign_rhs
    , regime     : reg_rhs
    , exponent   : exp_rhs
    , significand : significand_rhs})
  = rhs.sign_regime_exponent_significand() else { return rhs }


  ; let mult = Self::repr_mul_carry(significand, significand_rhs)
  ; let carry = 0 == (1 << Self::PRECISION-2) & mult
  ; let significand = mult << 1-carry as u32
  ; let power 
  = (regime as i16 + reg_rhs as i16) * Self::REGIME_BUMP as i16 
  + (exponent.0 + exp_rhs.0) as i16
  + carry as i16
  
  ; let mut regime =                  power / Self::REGIME_BUMP as i16
  ; let mut exponent = PositExponent((power % Self::REGIME_BUMP as i16) as i8)
  ; if exponent.0 < 0 { regime -= 1 ; exponent.0 += Self::REGIME_BUMP as i8 }


  ; let sign_bit = sign_bit ^ sign_rhs

  ; if  regime > (Self::PRECISION-2) as i16 { return if sign_bit{Self::MIN_NEG} else {Self::MAX_POS}; }
  ; if -regime > (Self::PRECISION-2) as i16 { return if sign_bit{Self::MAX_NEG} else {Self::MIN_POS}; }

  ; PositUnpacked::<<Self as PositRepr>::IRepr> { sign_bit, regime : regime as i8, exponent, significand}.encode() 
  }

  /// standard : "division(posit1,posit2) returns posit1÷posit2, rounded."
  pub const fn division(self : Self, rhs : Self) -> Self { self.multiplication(rhs.reciprocal()) }


  /// `NaR | 0 => NaR`, else `posit => 1/posit`
  ///
  /// non-standard
  pub const fn reciprocal(self : Self) -> Self 
  { if self.0 == Self::NAR.0 { return Self::NAR }
    Self(self.0.wrapping_neg() ^ (1<<Self::PRECISION-1))
  }

  /// `(posit, NaR | 0) | (NaR, posit) => None,` else `(posit, posit) => Some(Quot_Rem{quot : posit, rem : posit})`
  ///
  /// uses Donanld Knuth floor division for the quotient. `quot = floor(divisor / dividend), rem = divisor - divdend*quot`
  ///
  /// non-standard
  pub const fn quotient_remainder(self : Self, rhs : Self) -> Option<QuotRem<<Self as PositRepr>::IRepr>> 
  { if (self.0 == Self::NAR.0) | (rhs.0 == Self::NAR.0) | (rhs.0 == Self::ZERO.0) { return None }
  ; let quot = self.division(rhs).floor()
  ; let rem  = self.subtraction(rhs.multiplication(quot))
  ; Some(QuotRem{ quot, rem })
  }

  /// Decodes a posit into it's constituent parts
  const fn sign_regime_exponent_significand(self : Self) -> PositDecoded<<Self as PositRepr>::IRepr>
  { if self.0 == Self::NAR.0  { return PositDecoded::NaR }
    if self.0 == Self::ZERO.0 { return PositDecoded::Zero}

  ; let Posit(mut inner) = self
  ; let sign_bit = inner.is_negative()
  ; if sign_bit{inner = -inner}
  ; inner <<= 1

  ; let (regime, shift) = 
    if inner.is_negative() // leading 1 bit
    { let len = inner.leading_ones()  ; (len as i8-1, len) }
    else
    { let len = inner.leading_zeros() ; (-(len as i8) , len) }; // posit.abs() < 1
  ; inner <<= shift
  ; inner <<= 1

  ; let exponent = PositExponent((0b_11 & (inner >> Self::PRECISION-2)) as i8)

  ; let significand = (inner << 1  |  1 << Self::PRECISION-1) as  <Self as PositRepr>::URepr

  ; PositDecoded::NonZero(PositUnpacked { sign_bit, regime, exponent, significand})
  }
}


macro_rules! posit_methods
{($($TYPE:ty)*) => {$(
  impl Posit<$TYPE>
  { ////
    // 5.2 Basic functions of one posit value argument 
    ////
  
    /// +x => -x, -x => +x, else unchanged
    ///
    /// standard : "negate(posit) __ returns −posit.(7)"
    ///
    /// footnote : "(7) This is the 2’s complement of the posit representation. 2’s complement affects neither 0 nor NaR, since they are unsigned."
    pub const fn negate(self : Self) -> Self 
    {  Self( self.0.wrapping_neg() ) }
  
    /// -x => +x, else unchanged
    ///
    /// standard : "abs(posit) __ returns negate(posit) if posit < 0, else posit"
    pub const fn abs(self : Self) -> Self { Self( self.0.wrapping_abs() ) }
  
    /// +x => +1 , -x => -1, else unchanged
    ///
    /// standard : "sign(posit) __ returns a posit value: 1 if posit > 0, −1 if posit < 0, or 0 if posit = 0."
  
    pub const fn sign(self : Self) -> Self 
    { let nar_mask = !( ((self.0 == Self::NAR.0) as <Self as PositRepr>::IRepr) << Self::PRECISION - 2 ) // NAR => 0b_101..1 , _ => 0b_1..1
    ; Self( nar_mask & (self.0.signum() << Self::PRECISION - 2)) 
    }

    /// standard : "nearestInt(posit) __ returns the integer-valued posit value nearest to posit, and returns the nearest even integer-valued posit value if two integers are equally near."
    pub const fn nearest_int(self : Self) -> Self
    { if self.0 <= Self::FRACTION_LIMIT.negate().0 || self.0 >= Self::FRACTION_LIMIT.0 { return self } // note this includes NaR
      if self.0 >  Self::THREE_HALVES.negate().0   && self.0 <  Self::HALF.negate().0  { return Self::NEG_ONE }
      if self.0 >= Self::HALF.negate().0           && self.0 <= Self::HALF.0           { return Self::ZERO }
      if self.0 >  Self::HALF.0                    && self.0 <  Self::THREE_HALVES.0   { return Self::ONE }
      
      match self.sign_regime_exponent_significand() 
      { PositDecoded::NonZero(PositUnpacked { sign_bit, regime, exponent, significand })
        =>
        { let power = regime as u32 * Self::REGIME_BUMP + exponent.0 as u32 // this could be 0
        ; core::debug_assert!{ power <= Self::PRECISION-5 }
  
        ; let mask = !0 >> power+2
        ; let half_bit          = 0 != significand  &  mask+1
        ; let one_bit           = 0 != significand  &  mask+1 << 1
        ; let trailing_fraction = 0 != significand  &  mask
  
        ; let abs = self.abs()
  
        ; let integer = if half_bit & (!one_bit&trailing_fraction | one_bit)
          { abs.ceil() } else { abs.floor() }
  
        ; if sign_bit 
          { integer.negate() } else { integer }
        }
        PositDecoded::NaR | PositDecoded::Zero => { self } // this should in practice be unreachable
      }
  
    }

    /// int => int , +int.fract => +int + 1, -int.fract => -int, NAR => NAR
    ///
    /// standard : "ceil(posit) __ returns the smallest integer-valued posit value greater than or equal to posit."
    pub const fn ceil(self : Self) -> Self
    { if self.0 >  Self::NEG_ONE.0 && self.0 <= Self::ZERO.0 { return Self::ZERO }
      if self.0 >  Self::ZERO.0    && self.0 <= Self::ONE.0  { return Self::ONE  }
      if self.0 <= Self::FRACTION_LIMIT.negate().0 || self.0 >= Self::FRACTION_LIMIT.0 { return self } // note this includes NaR
      
      match self.sign_regime_exponent_significand()
      { PositDecoded::NonZero(PositUnpacked {sign_bit, regime, exponent, significand}) // at this point, posit.abs() > 1, so regime is  positive
        => 
        { let mut power = regime as u32 * Self::REGIME_BUMP + exponent.0 as u32
        ; let mask  = !(0 as <Self as PositRepr>::URepr) >> power+1
        ; let trailing_fraction = mask & significand
        
        ; let significand = 
          if sign_bit { significand & !mask }
          else if let Some(sig) = significand.checked_add(mask) { sig &  !mask }
          else { power+=1 ; 1 << Self::PRECISION-1 }
        
        ; let regime = (power / Self::REGIME_BUMP) as i8
        ; let exponent = PositExponent((power%Self::REGIME_BUMP) as i8)
        
        ; PositUnpacked::<<Self as PositRepr>::IRepr> {sign_bit, regime, exponent, significand}.encode()
        }
        PositDecoded::NaR | PositDecoded::Zero => { self } // this should in practice be unreachable
      }
    }

  
    /// int => int , +int.fract => int, -int.fract => -int - 1, NAR => NAR
    ///
    /// standard : "floor(posit) __ returns the largest integer-valued posit value less than or equal to posit."
    pub const fn floor(self : Self) -> Self
    { self.negate().ceil().negate() }
  
    /// The inverse of `Posit::prior`
    ///
    /// standard : "next(posit) __ returns the posit value of the lexicographic successor of posit’s representation. (8)"
    ///
    /// footnote : "(8) wrapping around, if necessary"
    ///
    /// note! lexicographical order is based on the total order; with wrapping these special cases should be in mind :
    ///
    /// MAX_POS => NAR, NAR => MIN_NEG, MAX_NEG => 0, 0 => MIN_POS
    pub const fn next(self : Self) -> Self{ Self( self.0.wrapping_add(1) ) }
  
    /// The inverse of `Posit::next`
    /// 
    /// standard : "prior(posit) __ returns the posit value of the lexicographic predecessor of posit’s representation. (9)"
    /// - footnote : "(9) wrapping around, if necessary"
    ///
    /// note! lexicographical order is based on the total order; with wrapping these special cases should be in mind : 
    ///
    /// MIN_NEG => NAR, NAR => MAX_POS, MIN_POS => 0, 0 => MAX_NEG
    pub const fn prior(self : Self) -> Self { Self( self.0.wrapping_sub(1) ) }

    ////
    // 5.3 Comparison functions of two posit value arguments
    ////
  
    pub const fn compare_equal         (self : Self, rhs : Self) -> bool {/*! standard : "compareEqual(posit1,posit2)        __ returns True if posit1 = posit2, else False." */ self.0 == rhs.0 }
    pub const fn compare_not_equal     (self : Self, rhs : Self) -> bool {/*! standard : "compareNotEqual(posit1,posit2)     __ returns True if posit1 ≠ posit2, else False." */ self.0 != rhs.0 }
    pub const fn compare_greater       (self : Self, rhs : Self) -> bool {/*! standard : "compareGreater(posit1,posit2)      __ returns True if posit1 > posit2, else False." */ self.0 >  rhs.0 }
    pub const fn compare_greater_equal (self : Self, rhs : Self) -> bool {/*! standard : "compareGreaterEqual(posit1,posit2) __ returns True if posit1 ≥ posit2, else False." */ self.0 >= rhs.0 }
    pub const fn compare_less          (self : Self, rhs : Self) -> bool {/*! standard : "compareLess(posit1,posit2)         __ returns True if posit1 < posit2, else False." */ self.0 <  rhs.0 }
    pub const fn compare_less_equal    (self : Self, rhs : Self) -> bool {/*! standard : "compareLessEqual(posit1,posit2)    __ returns True if posit1 ≤ posit2, else False." */ self.0 <= rhs.0 }

    ////
    // 5.4 Arithmetic functions of two posit value arguments
    ////

    /// standard : "addition(posit1,posit2) __ returns posit1+posit2, rounded."
    pub const fn addition(self : Self, rhs : Self) -> Self 
    { if self.0 == Self::NAR.0
      || rhs.0  == Self::NAR.0    { return Self::NAR }
    ; if self.0 == Self::ZERO.0   { return rhs}
    ; if rhs.0  == Self::ZERO.0   { return self}
    ; if self.0 == rhs.negate().0 { return Self::ZERO }
  
    let (big, small) = if self.0.abs() > rhs.0.abs() 
         {(self, rhs)} 
    else {(rhs, self)}    
  
    ; let PositDecoded::NonZero(PositUnpacked 
    { sign_bit   : sign_b
    , regime     : reg_b
    , exponent   : exp_b
    , significand : significand_b}) 
    = big.sign_regime_exponent_significand() else { core::unreachable!() }
    
    ; let PositDecoded::NonZero(PositUnpacked { 
        sign_bit   : sign_s
      , regime     : reg_s
      , exponent   : exp_s
      , significand : significand_s
      })
    = small.sign_regime_exponent_significand() else {  core::unreachable!() }
  
   
    // big small
    //  b+ s  =>   b+s
    // -b+-s  => -(b+s)
    //  b+-s  =>   b-s
    // -b+ s  => -(b-s)
  
    ; let power_b   = (reg_b as i16)*(Self::REGIME_BUMP as i16) + (exp_b.0 as i16)
    ; let power_s   = (reg_s as i16)*(Self::REGIME_BUMP as i16) + (exp_s.0 as i16)
  
    ; let power_diff = power_b - power_s
    ; if power_diff >= Self::PRECISION as i16 { return big }
  
    ; let significand_b_ = significand_b >> 1
    ; let significand_s_ = significand_s >> 1 >> power_diff
    
    ; let mut significand = if sign_b^sign_s
           { significand_b_ - significand_s_ }
      else { significand_b_ + significand_s_ }
    ; let leading = significand.leading_zeros()
    ; let power = power_b + 1 - leading as i16
    ; significand <<= leading

    ; let mut regime   =               (power / Self::REGIME_BUMP as i16) as i8
    ; let mut exponent = PositExponent((power % Self::REGIME_BUMP as i16) as i8)
    ; if exponent.0 < 0 { regime -=1 ; exponent.0 += Self::REGIME_BUMP as i8 }

    ; let sign_bit = sign_b
  
    ; PositUnpacked::<<Self as PositRepr>::IRepr> { sign_bit, regime, exponent, significand}.encode() 
    }
    /// standard : "subtraction(posit1,posit2) __ returns posit1−posit2, rounded."
    pub const fn subtraction(self : Self, rhs : Self) -> Self { self.addition(rhs.negate()) }
    /// standard : "multiplication(posit1,posit2) returns posit1×posit2, rounded."
    pub const fn multiplication(self : Self, rhs : Self) -> Self 
    { if rhs.0  == Self::NAR.0 { return Self::NAR }
    ; if rhs.0  == Self::ONE.0 { return self}
    ; if self.0 == Self::ONE.0 { return rhs}
    ; let PositDecoded::NonZero(PositUnpacked {sign_bit, regime, exponent, significand}) 
    = self.sign_regime_exponent_significand() else { return self }
  
    ; let PositDecoded::NonZero
      (PositUnpacked { 
        sign_bit   : sign_rhs
      , regime     : reg_rhs
      , exponent   : exp_rhs
      , significand : significand_rhs
      })
    = rhs.sign_regime_exponent_significand() else { return rhs }
  
  
    ; let mult = Self::repr_mul_carry(significand, significand_rhs)
    ; let carry = 0 == (1 << Self::PRECISION-2) & mult
    ; let significand = mult << 1-carry as u32
    ; let power 
    = (regime as i16 + reg_rhs as i16) * Self::REGIME_BUMP as i16 
    + (exponent.0 + exp_rhs.0) as i16
    + carry as i16
    
    ; let mut regime =                  power / Self::REGIME_BUMP as i16
    ; let mut exponent = PositExponent((power % Self::REGIME_BUMP as i16) as i8)
    ; if exponent.0 < 0 { regime -= 1 ; exponent.0 += Self::REGIME_BUMP as i8 }
  
  
    ; let sign_bit = sign_bit ^ sign_rhs
  
    ; if  regime > (Self::PRECISION-2) as i16 { return if sign_bit{Self::MIN_NEG} else {Self::MAX_POS}; }
    ; if -regime > (Self::PRECISION-2) as i16 { return if sign_bit{Self::MAX_NEG} else {Self::MIN_POS}; }
  
    ; PositUnpacked::<<Self as PositRepr>::IRepr> { sign_bit, regime : regime as i8, exponent, significand}.encode() 
    }

    /// standard : "division(posit1,posit2) returns posit1÷posit2, rounded."
    pub const fn division(self : Self, rhs : Self) -> Self { self.multiplication(rhs.reciprocal()) }

    /// `NaR | 0 => NaR`, else `posit => 1/posit`
    ///
    /// non-standard
    pub const fn reciprocal(self : Self) -> Self 
    { if self.0 == Self::NAR.0 { return Self::NAR }
      Self(self.0.wrapping_neg() ^ (1<<Self::PRECISION-1))
    }
  
    /// `(posit, NaR | 0) | (NaR, posit) => None,` else `(posit, posit) => Some(Quot_Rem{quot : posit, rem : posit})`
    ///
    /// uses Donanld Knuth floor division for the quotient. `quot = floor(divisor / dividend), rem = divisor - divdend*quot`
    ///
    /// non-standard
    pub const fn quotient_remainder(self : Self, rhs : Self) -> Option<QuotRem<<Self as PositRepr>::IRepr>> 
    { if (self.0 == Self::NAR.0) | (rhs.0 == Self::NAR.0) | (rhs.0 == Self::ZERO.0) { return None }
    ; let quot = self.division(rhs).floor()
    ; let rem  = self.subtraction(rhs.multiplication(quot))
    ; Some(QuotRem{ quot, rem })
    }

    const fn sign_regime_exponent_significand(self : Self) -> PositDecoded<<Self as PositRepr>::IRepr>
    { if self.0 == Self::NAR.0  { return PositDecoded::NaR }
      if self.0 == Self::ZERO.0 { return PositDecoded::Zero}
  
    ; let Posit(mut inner) = self
    ; let sign_bit = inner.is_negative()
    ; if sign_bit{inner = -inner}
    ; inner <<= 1
  
    ; let (regime, shift) = 
      if inner.is_negative() // leading 1 bit
           { let len = inner.leading_ones()  ; (  len as i8-1, len) }
      else { let len = inner.leading_zeros() ; (-(len as i8) , len) } // posit.abs() < 1
    ; inner <<= shift
    ; inner <<= 1
  
    ; let exponent = PositExponent((0b_11 & (inner >> Self::PRECISION-2)) as i8)
  
    ; let significand = (inner << 1  |  1 << Self::PRECISION-1) as  <Self as PositRepr>::URepr
  
    ; PositDecoded::NonZero(PositUnpacked { sign_bit, regime, exponent, significand})
    }
  }
)*};}
posit_methods!{i8 i16 i32 isize i128}



macro_rules! tmp_posit_dbg_fmt 
{($($TYPE:ty)*) => {$(
impl  core::fmt::Debug for Posit<$TYPE>
{ fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result 
  { // core::write!( f, "{:#x}:{:?}", self.0,self.sign_regime_exponent_significand() )
  ; let name = core::any::type_name::<$TYPE>()
  ; match self.sign_regime_exponent_significand() 
    { PositDecoded::NonZero(PositUnpacked{ sign_bit, regime, exponent, significand } )
                         => { let as_float = if sign_bit {-1.0}else{1.0}*(significand as f64 * 2.0f64.powf(4.0*regime as f64 + exponent.0 as f64 - (<Self as PositRepr>::URepr::BITS-1) as f64) )
                            ; core::write!( f, "[Posit<{}> {:#x} [{} {} {} {:#x}] {}_f64]",name, self.0,if sign_bit {'-'} else {'+'}, regime, exponent.0, significand, as_float)}
      PositDecoded::Zero => { core::write!( f, "[Posit<{}> Zero]", name ) }
      PositDecoded::NaR  => { core::write!( f, "[Posit<{}> NaR]", name) } 
    }
  }}  
)*};}
tmp_posit_dbg_fmt!{i8 i16 i32 i64 i128 isize}

////
// Format specific impl
///

// this should only be templated for non i128 types
impl Posit<i64>
{ /// used internally for multiplying
  const fn repr_mul_carry(left : <Self as PositRepr>::URepr, right :  <Self as PositRepr>::URepr)-> <Self as PositRepr>::URepr
  { ( left  as <Self as PositRepr>::UWide 
    * right as <Self as PositRepr>::UWide
    >> Self::PRECISION
    ) as <Self as PositRepr>::URepr
  }
}

macro_rules! posit_non_wide
{($($TYPE:ty)*) => {$(
  impl Posit<$TYPE>
  { /// used internally for multiplying
    const fn repr_mul_carry(left : <Self as PositRepr>::URepr, right :  <Self as PositRepr>::URepr)-> <Self as PositRepr>::URepr
    { ( left  as <Self as PositRepr>::UWide 
      * right as <Self as PositRepr>::UWide 
      >> Self::PRECISION
      ) as <Self as PositRepr>::URepr
    }
  }
)*};}
posit_non_wide!{i8 i16 i32}


impl Posit<i128>
{ /// a specific implementation used internally for multiplying the fraction, made just for i128 as it does not have a unsigned wide associated type
  const fn repr_mul_carry(left : <Self as PositRepr>::URepr, right :  <Self as PositRepr>::URepr)-> <Self as PositRepr>::URepr
  { let half_bits = Self::PRECISION / 2
  ; let bot_mask = !(!0 << Self::PRECISION / 2)
  ; let [lt, lb] = [left  >> half_bits, left  & bot_mask]
  ; let [rt, rb] = [right >> half_bits, right & bot_mask]

  ; let [lt_x_rt, lt_x_rb, lb_x_rt, lb_x_rb] = [ lt*rt, lt*rb, lb*rt, lb*rb ]

  ; let low_bits_sum = (lb_x_rb >> half_bits) + (lt_x_rb & bot_mask) + (lb_x_rt & bot_mask)
  ; lt_x_rt + (lt_x_rb >> half_bits) + (lb_x_rt >> half_bits) + (low_bits_sum >> half_bits)
  }
}
macro_rules! posit_wide_impl
{($($TYPE:ty)*) => {$(

  impl Posit<$TYPE>
  { /// a specific implementation used internally for multiplying the fraction, made just for i128 as it does not have a unsigned wide associated type
    const fn repr_mul_carry(left : <Self as PositRepr>::URepr, right :  <Self as PositRepr>::URepr)-> <Self as PositRepr>::URepr
    { let half_bits = Self::PRECISION / 2
    ; let bot_mask = !(!0 << Self::PRECISION / 2)
    ; let [lt, lb] = [left  >> half_bits, left  & bot_mask]
    ; let [rt, rb] = [right >> half_bits, right & bot_mask]
  
    ; let [lt_x_rt, lt_x_rb, lb_x_rt, lb_x_rb] = [ lt*rt, lt*rb, lb*rt, lb*rb ]
  
    ; let low_bits_sum = (lb_x_rb >> half_bits) + (lt_x_rb & bot_mask) + (lb_x_rt & bot_mask)
    ; lt_x_rt + (lt_x_rb >> half_bits) + (lb_x_rt >> half_bits) + (low_bits_sum >> half_bits)
    }
  }
)*};}
posit_wide_impl!{isize}

#[test]#[cfg(test)] 
fn test_repr_mul_carry_isize()
{ type P = Posit<i64>
; let (l, r) = ((1 << P::PRECISION-1) | 1, 1 << P::PRECISION-1); let expected = 1 << P::PRECISION-2 // note the trailing `| 1` gets cut off
// ; extern crate std; std::println!("left    : {l:#0x}\nright   : {r:#0x}\nexpected: {expected:#0x}")
; core::assert_eq!{P::repr_mul_carry(l,r), expected}
}
macro_rules! posit_unit_tests 
{($($TYPE:ty)*) => {
  #[cfg(test)]#[test]
  fn test_repr_mul_carry_templated() 
  {$({
    type P = Posit<$TYPE>
    ; let (l, r) = ((1 << P::PRECISION-1) | 1, 1 << P::PRECISION-1)
    ; let expected = 1 << P::PRECISION-2 // note the trailing `| 1` gets cut off
    ; core::assert_eq!{P::repr_mul_carry(l,r), expected}
  }  )*}

  #[cfg(test)]#[test]
  fn test_sign_regime_exponent_significand()
  {$({ use PositUnpacked as PU
  ; type Repr = $TYPE
  ; type P = Posit<Repr>
  ; const P : fn(Repr)->P = Posit::<Repr>
  ; let nz = PositDecoded::<Repr>::NonZero
  ; let sres = P::sign_regime_exponent_significand
  ;
    core::assert_eq!{ sres(P::MAX_POS), nz(PU { sign_bit : false, regime : (Repr::BITS-2) as i8   , exponent : PositExponent(0), significand : 1 << Repr::BITS-1}) }
    core::assert_eq!{ sres(P::ONE    ), nz(PU { sign_bit : false, regime : 0                      , exponent : PositExponent(0), significand : 1 << Repr::BITS-1}) }
    core::assert_eq!{ sres(P::MIN_POS), nz(PU { sign_bit : false, regime : -((Repr::BITS-2) as i8), exponent : PositExponent(0), significand : 1 << Repr::BITS-1}) }
    core::assert_eq!{ sres(P::ZERO   ), PositDecoded::Zero }
    core::assert_eq!{ sres(P::MAX_NEG), nz(PU { sign_bit : true , regime : -((Repr::BITS-2) as i8), exponent : PositExponent(0), significand : 1 << Repr::BITS-1}) }
    core::assert_eq!{ sres(P::NEG_ONE), nz(PU { sign_bit : true , regime : 0                      , exponent : PositExponent(0), significand : 1 << Repr::BITS-1}) }
    core::assert_eq!{ sres(P::MIN_NEG), nz(PU { sign_bit : true , regime : (Repr::BITS-2) as i8   , exponent : PositExponent(0), significand : 1 << Repr::BITS-1}) }
    core::assert_eq!{ sres(P::NAR    ), PositDecoded::NaR }

    core::assert_eq!
    { sres(P(0b_010001 << Repr::BITS-6)) // 1 + 1/2  
    , nz(PU { sign_bit : false , regime : 0, exponent : PositExponent(0), significand : 0b_11 << Repr::BITS-2}) 
    }
    core::assert_eq!
    { sres(P(0b_010001 << Repr::BITS-6).negate()) // -(1 + 1/2)
    , nz(PU { sign_bit : true , regime : 0, exponent : PositExponent(0), significand : 0b_11 << Repr::BITS-2}) 
    }
    
    core::assert_eq!
    { sres(P(0b_00111 << Repr::BITS-5)) // 1/2  
    , nz(PU { sign_bit : false , regime : -1, exponent : PositExponent(0b_11), significand : 0b_1 << Repr::BITS-1}) 
    }
    core::assert_eq!
    { sres(P(0b_00111 << Repr::BITS-5).negate()) // - 1/2  
    , nz(PU { sign_bit : true , regime : -1, exponent : PositExponent(0b_11), significand : 0b_1 << Repr::BITS-1}) 
    }
    
  })*}

  #[cfg(test)]#[test] fn test_negate()
  {$({
    type P = Posit<$TYPE>;
    core::assert_eq!{ P::NAR    .negate() , P::NAR}
    core::assert_eq!{ P::MAX_POS.negate() , P::MIN_NEG}
    core::assert_eq!{ P::ONE    .negate() , P::NEG_ONE}
    core::assert_eq!{ P::MIN_POS.negate() , P::MAX_NEG}
    core::assert_eq!{ P::ZERO   .negate() , P::ZERO}
    core::assert_eq!{ P::MAX_NEG.negate() , P::MIN_POS}
    core::assert_eq!{ P::NEG_ONE.negate() , P::ONE}
    core::assert_eq!{ P::MIN_NEG.negate() , P::MAX_POS}
  })*}

  #[cfg(test)]#[test] fn test_abs()
  {$({
  type P = Posit<$TYPE>;
  core::assert_eq!{ P::NAR    .abs() , P::NAR}
  core::assert_eq!{ P::MAX_POS.abs() , P::MAX_POS}
  core::assert_eq!{ P::ONE    .abs() , P::ONE}
  core::assert_eq!{ P::MIN_POS.abs() , P::MIN_POS}
  core::assert_eq!{ P::ZERO   .abs() , P::ZERO}
  core::assert_eq!{ P::MAX_NEG.abs() , P::MIN_POS}
  core::assert_eq!{ P::NEG_ONE.abs() , P::ONE}
  core::assert_eq!{ P::MIN_NEG.abs() , P::MAX_POS}
  })*}
  

  #[cfg(test)]#[test] fn test_sign()
  {$({  
    type P = Posit<$TYPE>;
    core::assert_eq!{ P::NAR    .sign() , P::NAR}
    core::assert_eq!{ P::MAX_POS.sign() , P::ONE}
    core::assert_eq!{ P::ONE    .sign() , P::ONE}
    core::assert_eq!{ P::MIN_POS.sign() , P::ONE}
    core::assert_eq!{ P::ZERO   .sign() , P::ZERO}
    core::assert_eq!{ P::MAX_NEG.sign() , P::NEG_ONE}
    core::assert_eq!{ P::NEG_ONE.sign() , P::NEG_ONE}
    core::assert_eq!{ P::MIN_NEG.sign() , P::NEG_ONE}
  })*}


  #[cfg(test)]#[test] fn test_ceil()
  {$({  
    type Repr = $TYPE;
    type P = Posit<Repr>;
    const P : fn(Repr)->P = Posit::<Repr>::from_bits;
    core::assert_eq!{ P::NAR    .ceil(), P::NAR}
    core::assert_eq!{ P::MAX_POS.ceil(), P::MAX_POS}
    core::assert_eq!{ P::ONE    .ceil(), P::ONE}
    core::assert_eq!{ P::MIN_POS.ceil(), P::ONE}
    core::assert_eq!{ P::ZERO   .ceil(), P::ZERO}
    core::assert_eq!{ P::MAX_NEG.ceil(), P::ZERO}
    core::assert_eq!{ P::NEG_ONE.ceil(), P::NEG_ONE}
    core::assert_eq!{ P::MIN_NEG.ceil(), P::MIN_NEG}
    
    core::assert_eq!{ P::HALF           .ceil(), P::ONE} 
    core::assert_eq!{ P::HALF.negate()  .ceil(), P::ZERO}
    core::assert_eq!{ P::ONE.next()     .ceil(), P(0b_01001 << Repr::BITS-5)}
    core::assert_eq!{ P::NEG_ONE.prior().ceil(), P::NEG_ONE}                 
    
    core::assert_eq!{ P(0b_0100111 << Repr::BITS-7).ceil(), P(0b_0101 << Repr::BITS-4) }
    if Repr::BITS > 8 {core::assert_eq!{ P((0b_01011_1111_i128 << Repr::BITS.wrapping_sub(9)) as Repr).ceil(), P(0b_011 << Repr::BITS-3) }};
  })*}


  #[cfg(test)]#[test] fn test_floor()
  {$({  
    type Repr = $TYPE;
    type P = Posit<Repr>;
    const P : fn(Repr)->P = Posit::<Repr>::from_bits;
    core::assert_eq!{ P::NAR    .floor(), P::NAR}
    core::assert_eq!{ P::MAX_POS.floor(), P::MAX_POS}
    core::assert_eq!{ P::ONE    .floor(), P::ONE}
    core::assert_eq!{ P::MIN_POS.floor(), P::ZERO}
    core::assert_eq!{ P::ZERO   .floor(), P::ZERO}
    core::assert_eq!{ P::MAX_NEG.floor(), P::NEG_ONE}
    core::assert_eq!{ P::NEG_ONE.floor(), P::NEG_ONE}
    core::assert_eq!{ P::MIN_NEG.floor(), P::MIN_NEG}
    
    core::assert_eq!{ P::HALF           .floor(), P::ZERO}
    core::assert_eq!{ P::HALF.negate()  .floor(), P::NEG_ONE}
    core::assert_eq!{ P::ONE.next()     .floor(), P::ONE}
    core::assert_eq!{ P::NEG_ONE.prior().floor(), P(0b_01001 << Repr::BITS-5).negate()}
    
    core::assert_eq!{ P(0b_0100111 << Repr::BITS-7).floor(), P(0b_010011 << Repr::BITS-6) }
    if Repr::BITS > 8 {core::assert_eq!{ 
      P((0b_01011_1111_i128 << Repr::BITS.wrapping_sub(9)) as Repr).floor()
    , P((0b_01011_111_i128  << Repr::BITS-8) as Repr) 
    }}
  ;
  })*}

  #[cfg(test)]#[test] fn test_next()
  {$({    
    type P = Posit<$TYPE>;
    core::assert_eq!{ P::MAX_POS.next() , P::NAR}
    core::assert_eq!{ P::NAR    .next() , P::MIN_NEG}
    core::assert_eq!{ P::ZERO   .next() , P::MIN_POS}
    core::assert_eq!{ P::MAX_NEG.next() , P::ZERO}
  })*}

  #[cfg(test)]#[test] fn test_prior()
  {$({    
    type P = Posit<$TYPE>;
    core::assert_eq!{ P::MIN_NEG.prior() , P::NAR}
    core::assert_eq!{ P::NAR    .prior() , P::MAX_POS}
    core::assert_eq!{ P::MIN_POS.prior() , P::ZERO}
    core::assert_eq!{ P::ZERO   .prior() , P::MAX_NEG}
  })*}

  #[cfg(test)]#[test] fn test_nearest_int()
  {$({
    type Repr = $TYPE;
    type P = Posit<Repr>;
    const P : fn(Repr)->P = Posit::<Repr>::from_bits;
    core::assert_eq!{ P::NAR    .nearest_int() , P::NAR}
    core::assert_eq!{ P::MAX_POS.nearest_int() , P::MAX_POS}
    core::assert_eq!{ P::ONE    .nearest_int() , P::ONE}
    core::assert_eq!{ P::MIN_POS.nearest_int() , P::ZERO}
    core::assert_eq!{ P::ZERO   .nearest_int() , P::ZERO}
    core::assert_eq!{ P::MAX_NEG.nearest_int() , P::ZERO}
    core::assert_eq!{ P::NEG_ONE.nearest_int() , P::NEG_ONE}
    core::assert_eq!{ P::MIN_NEG.nearest_int() , P::MIN_NEG}

    core::assert_eq!{ P::HALF                .nearest_int() , P::ZERO }
    core::assert_eq!{ P::HALF       .negate().nearest_int() , P::ZERO }
    core::assert_eq!{ P::HALF.next()         .nearest_int() , P::ONE }
    core::assert_eq!{ P::HALF.next().negate().nearest_int() , P::NEG_ONE }

    core::assert_eq!{ P::THREE_HALVES                 .nearest_int() , P::THREE_HALVES.ceil() }
    core::assert_eq!{ P::THREE_HALVES        .negate().nearest_int() , P::THREE_HALVES.ceil().negate() }
    core::assert_eq!{ P::THREE_HALVES.next()          .nearest_int() , P::THREE_HALVES.ceil() }
    core::assert_eq!{ P::THREE_HALVES.next() .negate().nearest_int() , P::THREE_HALVES.ceil().negate() }
    core::assert_eq!{ P::THREE_HALVES.prior()         .nearest_int() , P::THREE_HALVES.floor() }
    core::assert_eq!{ P::THREE_HALVES.prior().negate().nearest_int() , P::THREE_HALVES.floor().negate() }

    core::assert_eq!{ P::FRACTION_LIMIT.prior()         .nearest_int() , P::FRACTION_LIMIT } 
    core::assert_eq!{ P::FRACTION_LIMIT.prior().negate().nearest_int() , P::FRACTION_LIMIT.negate() }
    core::assert_ne!{ P::FRACTION_LIMIT.next()          .nearest_int() , P::FRACTION_LIMIT } 
    core::assert_ne!{ P::FRACTION_LIMIT.next() .negate().nearest_int() , P::FRACTION_LIMIT.negate() }
  })*}


  #[cfg(test)]#[test] fn test_reciprocal()
  {$({
    type Repr = $TYPE;
    type P = Posit<Repr>;
    const P : fn(Repr)->P = Posit::<Repr>::from_bits;
    core::assert_eq!{ P::NAR    .reciprocal() , P::NAR}
    core::assert_eq!{ P::MAX_POS.reciprocal() , P::MIN_POS}
    core::assert_eq!{ P::ONE    .reciprocal() , P::ONE}
    core::assert_eq!{ P::MIN_POS.reciprocal() , P::MAX_POS}
    core::assert_eq!{ P::ZERO   .reciprocal() , P::NAR}
    core::assert_eq!{ P::MAX_NEG.reciprocal() , P::MIN_NEG}
    core::assert_eq!{ P::NEG_ONE.reciprocal() , P::NEG_ONE}
    core::assert_eq!{ P::MIN_NEG.reciprocal() , P::MAX_NEG}

    core::assert_eq!{ P::HALF         .reciprocal() , P(0b_01001 << Repr::BITS-5) }
    core::assert_eq!{ P::HALF.negate().reciprocal() , P(0b_01001 << Repr::BITS-5).negate() }
    core::assert_eq!{ P(0b_01001 << Repr::BITS-5)         .reciprocal() , P::HALF }
    core::assert_eq!{ P(0b_01001 << Repr::BITS-5).negate().reciprocal() , P::HALF.negate() }

    core::assert_eq!{ P::ONE.next() .reciprocal() , P::ONE.prior() }
    core::assert_eq!{ P::ONE.prior().reciprocal() , P::ONE.next() }
    core::assert_eq!{ P::NEG_ONE.next() .reciprocal() , P::NEG_ONE.prior() }
    core::assert_eq!{ P::NEG_ONE.prior().reciprocal() , P::NEG_ONE.next() }

   
  ; const RANGE : core::ops::RangeInclusive<Repr> = Repr::MIN..=Repr::MAX 
  ; let mut misses = 0;
    ; if Repr::BITS <= 16 
      { for each in RANGE 
        { if let P::ZERO | P::NAR = P(each) { continue }
        ; core::assert_eq!{ P(each).reciprocal().reciprocal(), P(each) }
        }
      }
  })*}

  #[cfg(test)]#[test] fn test_multiplication()
  {$({
    type Repr = $TYPE
  ; type P = Posit<Repr>
  ; const P : fn(Repr)->P = Posit::<Repr>::from_bits
  ; let edge_cases = [P::NAR, P::MIN_NEG, P::NEG_ONE, P::MAX_NEG, P::ZERO, P::MIN_POS, P::ONE, P::MAX_POS]

  ; for idx in 0..edge_cases.len() 
    { core::assert_eq!{ P::NAR.multiplication(edge_cases[idx]), P::NAR} 
      core::assert_eq!{ edge_cases[idx].multiplication(P::NAR), P::NAR}
     
     
      core::assert_eq!{ P::ONE.multiplication(edge_cases[idx]), edge_cases[idx]}
      core::assert_eq!{ edge_cases[idx].multiplication(P::ONE), edge_cases[idx]}
     
      core::assert_eq!{ P::NEG_ONE.multiplication(edge_cases[idx]), edge_cases[idx].negate()}
      core::assert_eq!{ edge_cases[idx].multiplication(P::NEG_ONE), edge_cases[idx].negate()}
       
    ; if edge_cases[idx] != P::NAR
      { 
      core::assert_eq!{ P::ZERO.multiplication(edge_cases[idx]), P::ZERO}
      core::assert_eq!{ edge_cases[idx].multiplication(P::ZERO), P::ZERO}
      }
    } 

  ; const RANGE : core::ops::RangeInclusive<Repr> = Repr::MIN..=Repr::MAX 
  ; if Repr::BITS == 8 
    { for left in RANGE {for right in RANGE {  
      core::assert_eq!{ P(left).multiplication(P(right)), P(right).multiplication(P(left)) }
    }}}


  ; if Repr::BITS <= 16
    { for each in RANGE 
      { if let P::ZERO | P::NAR = P(each) { continue }
      ; let mul_by_recip = P(each).multiplication(P(each).reciprocal())
      ; core::assert!
        {  P::ONE.0       <= mul_by_recip.0 
        && mul_by_recip.0 <= P(0b_0100_0001<<Repr::BITS-8).0 /* Posit(9/8) */ 
        }
      }
    }

  ; let two = P(0b_01001<<Repr::BITS-5)
  ; core::assert_eq!{ two   .multiplication(P::ONE), two }
    core::assert_eq!{ P::ONE.multiplication(two)   , two }

  ; let four = P(0b_0101<<Repr::BITS-4)
  ; core::assert_eq!{ two.multiplication(two)      , four }
    core::assert_eq!{ four.multiplication(four)    , P(0b_011<<Repr::BITS-3)}
    
  ; let fourth = four.reciprocal()
  ; core::assert_eq!{ fourth.multiplication(fourth.negate()).reciprocal(), four.multiplication(four).negate() }
  
  })*}

  #[cfg(test)]#[test] fn test_addition()
  {$({
    type Repr = $TYPE
  ; type P = Posit<Repr>
  ; const P : fn(Repr)->P = Posit::<Repr>::from_bits
  ; let edge_cases = [P::NAR, P::MIN_NEG, P::NEG_ONE, P::MAX_NEG, P::ZERO, P::MIN_POS, P::ONE, P::MAX_POS]
  
  ; for idx in 0..edge_cases.len() 
    { core::assert_eq!{ P::NAR.addition(edge_cases[idx]), P::NAR}
      core::assert_eq!{ edge_cases[idx].addition(P::NAR), P::NAR}

      core::assert_eq!{ P::ZERO.addition(edge_cases[idx]), edge_cases[idx]}
      core::assert_eq!{ edge_cases[idx].addition(P::ZERO), edge_cases[idx]}
    }
  
  ; const RANGE : core::ops::RangeInclusive<Repr> = Repr::MIN..=Repr::MAX 
  ; if Repr::BITS == 8 
    { for left in RANGE {for right in RANGE {  
      core::assert_eq!{ P(left).addition(P(right)), P(right).addition(P(left)) }
    }}}
  ; if Repr::BITS <= 16
    { for each in RANGE.into_iter().map(P)
      { if each.0 == P::NAR.0 { continue }
      ; core::assert_eq!{ each.addition(each.negate()), P::ZERO }

      ; let two          = P(0b_01001<<Repr::BITS-5)
      ; let regime_power = P(0b_0110<<Repr::BITS-4) // Posit(2^4)
      ; let limit        = P::FRACTION_LIMIT.division(regime_power.multiplication(regime_power)) // regime.pow(es) ?

      ; if   each.0 >  limit.negate().0
          && each.0 <= limit.negate().reciprocal().0
        ||   each.0 == P::ZERO.0
        ||   each.0 >= limit.reciprocal().0
          && each.0 <  limit.0
        { core::assert_eq!{ each.addition(each), each.multiplication(two) }}
      }
    }

  ; core::assert_eq!{ P::HALF.addition(P::ONE), P::THREE_HALVES }
  ; core::assert_eq!{ P::ONE.addition(P::HALF.negate()), P::HALF }

  // 343.5 + 57.75 => 401.25
  ; if Repr::BITS > 16 
    { core::assert_eq!
      {           P(((0b_0_1110_00__01010111_1  % Repr::MAX as i128)as Repr)<<Repr::BITS-16)
        .addition(P(((0b_0_110_01__11001_11     % Repr::MAX as i128)as Repr)<<Repr::BITS-13))
      , P(((0b_0_1110_00__10010001_01 % Repr::MAX as i128)as Repr)<<Repr::BITS-17)
      }
    }
  ; if Repr::BITS == 8
    { core::assert_eq!{P::ONE.next().addition(P::ONE.negate()), P(0b_0_01_01<<Repr::BITS-5)}
    }
  ; core::assert_eq!{ P::ONE.addition(P::ONE), P(0b_01001<<Repr::BITS-5) }
  })*}

  #[cfg(test)]#[test] fn test_quotient_remainder()
  {$({ 
    type Repr = $TYPE
  ; type P = Posit<Repr>
  ; const P : fn(Repr)->P = Posit::<Repr>::from_bits
  ; let edge_cases = [P::NAR, P::MIN_NEG, P::NEG_ONE, P::MAX_NEG, P::ZERO, P::MIN_POS, P::ONE, P::MAX_POS]
  ; for idx in 0..edge_cases.len() 
    { core::assert_eq!{ P::NAR.quotient_remainder(edge_cases[idx]), None}
      core::assert_eq!{ edge_cases[idx].quotient_remainder(P::NAR), None}
      core::assert_eq!{ edge_cases[idx].quotient_remainder(P::ZERO), None}
    }

    core::assert_eq!{ P::MAX_POS.quotient_remainder(P::MIN_POS), Some(QuotRem{quot : P::MAX_POS, rem : P::MAX_POS}) }
    core::assert_eq!{ P::MAX_POS.quotient_remainder(P::MAX_POS), Some(QuotRem{quot : P::ONE    , rem : P::ZERO}) }
    core::assert_eq!{ P::MIN_POS.quotient_remainder(P::MAX_POS), Some(QuotRem{quot : P::ZERO   , rem : P::MIN_POS}) }
    core::assert_eq!{ P::MIN_POS.quotient_remainder(P::MIN_POS), Some(QuotRem{quot : P::ONE    , rem : P::ZERO}) }

    core::assert_eq!{ ( P::ONE + P::ONE + P::ONE).quotient_remainder(P::ONE + P::ONE), Some(QuotRem{quot : P::ONE             , rem : P::ONE}) }
    core::assert_eq!{ ( P::THREE_HALVES)         .quotient_remainder(P::ONE)         , Some(QuotRem{quot : P::ONE             , rem : P::HALF}) }
    core::assert_eq!{ (-P::THREE_HALVES)         .quotient_remainder(P::ONE)         , Some(QuotRem{quot : -(P::ONE + P::ONE) , rem : P::HALF}) }
    core::assert_eq!{ ( P::THREE_HALVES)         .quotient_remainder(P::NEG_ONE)     , Some(QuotRem{quot : -(P::ONE + P::ONE) , rem : -P::HALF}) }
    
  })*}


  #[cfg(test)]#[test] fn test_dbg_posit()
  {{ 
    type Repr = i16
  ; type P = Posit<Repr>
  ; const P : fn(Repr)->P = Posit::<Repr>::from_bits
  ; for each in Repr::MIN..=Repr::MAX
    { extern crate std
    ; std::println!("{:?}", P(each))
    }
  }}
};}

posit_unit_tests!{i8 i16 i32 i64 i128 isize}

