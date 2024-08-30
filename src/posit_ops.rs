//! Rust operator trait implementations for Posit<T>
use super::Posit;
use crate::prelude::*;
use core::ops::
{ Neg
, Add, AddAssign
, Mul, MulAssign
, Sub, SubAssign
, Div, DivAssign
, Rem, RemAssign
};
use core::cmp::{PartialOrd, Ord};

macro_rules! impl_operator_traits {($($TYPE:ty)*) => {$(
  impl Ord for Posit<$TYPE> 
  { fn cmp(&self, other: &Self) -> core::cmp::Ordering { self.0.cmp(&other.0) } }
  
  impl PartialOrd for Posit<$TYPE> 
  { fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> { self.0.partial_cmp(&other.0)}}
  
  impl Neg for Posit<$TYPE>  { type Output = Self      ; fn neg(self) -> Self::Output { self.negate() } }
  impl Neg for &Posit<$TYPE> { type Output = Posit<$TYPE>; fn neg(self) -> Self::Output { self.negate() } }

  impl Add                for  Posit<$TYPE> { type Output = Self         ; fn add(self, rhs: Self         ) -> Self::Output { self.addition( rhs) } }
  impl Add                for &Posit<$TYPE> { type Output = Posit<$TYPE> ; fn add(self, rhs: Self         ) -> Self::Output { self.addition(*rhs) } }
  impl Add<Posit<$TYPE>>  for &Posit<$TYPE> { type Output = Posit<$TYPE> ; fn add(self, rhs: Posit<$TYPE> ) -> Self::Output { self.addition( rhs) } }
  impl Add<&Posit<$TYPE>> for  Posit<$TYPE> { type Output = Self         ; fn add(self, rhs: &Posit<$TYPE>) -> Self::Output { self.addition(*rhs) } }

  impl AddAssign                for Posit<$TYPE> { fn add_assign(&mut self, rhs: Self)          { *self = self.addition( rhs)  } }
  impl AddAssign<&Posit<$TYPE>> for Posit<$TYPE> { fn add_assign(&mut self, rhs: &Posit<$TYPE>) { *self = self.addition(*rhs)  } }

  impl Mul                for  Posit<$TYPE> { type Output = Self         ; fn mul(self, rhs: Self         ) -> Self::Output { self.multiplication( rhs) } }
  impl Mul                for &Posit<$TYPE> { type Output = Posit<$TYPE> ; fn mul(self, rhs: Self         ) -> Self::Output { self.multiplication(*rhs) } }
  impl Mul<Posit<$TYPE>>  for &Posit<$TYPE> { type Output = Posit<$TYPE> ; fn mul(self, rhs: Posit<$TYPE> ) -> Self::Output { self.multiplication( rhs) } }
  impl Mul<&Posit<$TYPE>> for  Posit<$TYPE> { type Output = Self         ; fn mul(self, rhs: &Posit<$TYPE>) -> Self::Output { self.multiplication(*rhs) } }
  
  impl MulAssign                for Posit<$TYPE> { fn mul_assign(&mut self, rhs: Self)          { *self = self.multiplication( rhs)  } }
  impl MulAssign<&Posit<$TYPE>> for Posit<$TYPE> { fn mul_assign(&mut self, rhs: &Posit<$TYPE>) { *self = self.multiplication(*rhs)  } }

  impl Sub                for  Posit<$TYPE> { type Output = Self         ; fn sub(self, rhs: Self          ) -> Self::Output { self.subtraction( rhs) } }
  impl Sub                for &Posit<$TYPE> { type Output = Posit<$TYPE> ; fn sub(self, rhs: Self          ) -> Self::Output { self.subtraction(*rhs) } }
  impl Sub<Posit<$TYPE>>  for &Posit<$TYPE> { type Output = Posit<$TYPE> ; fn sub(self, rhs: Posit<$TYPE>  ) -> Self::Output { self.subtraction( rhs) } }
  impl Sub<&Posit<$TYPE>> for  Posit<$TYPE> { type Output = Self         ; fn sub(self, rhs: &Posit<$TYPE> ) -> Self::Output { self.subtraction(*rhs) } }

  impl SubAssign                for Posit<$TYPE> { fn sub_assign(&mut self, rhs: Self)          { *self = self.subtraction( rhs)  } }
  impl SubAssign<&Posit<$TYPE>> for Posit<$TYPE> { fn sub_assign(&mut self, rhs: &Posit<$TYPE>) { *self = self.subtraction(*rhs)  } }

  impl Div                for  Posit<$TYPE> { type Output = Self         ; fn div(self, rhs: Self         ) -> Self::Output { self.division( rhs) } }
  impl Div                for &Posit<$TYPE> { type Output = Posit<$TYPE> ; fn div(self, rhs: Self         ) -> Self::Output { self.division(*rhs) } }
  impl Div<Posit<$TYPE>>  for &Posit<$TYPE> { type Output = Posit<$TYPE> ; fn div(self, rhs: Posit<$TYPE> ) -> Self::Output { self.division( rhs) } }
  impl Div<&Posit<$TYPE>> for  Posit<$TYPE> { type Output = Self         ; fn div(self, rhs: &Posit<$TYPE>) -> Self::Output { self.division(*rhs) } }

  impl DivAssign                for Posit<$TYPE> { fn div_assign(&mut self, rhs: Self)          { *self = self.division( rhs)  } }
  impl DivAssign<&Posit<$TYPE>> for Posit<$TYPE> { fn div_assign(&mut self, rhs: &Posit<$TYPE>) { *self = self.division(*rhs)  } }

  impl Rem                for  Posit<$TYPE> { type Output = Self         ; fn rem(self, rhs:  Self      ) -> Self::Output { self.quotient_remainder( rhs).map_or( Self::NAR            , |q_r| q_r.rem) } }
  impl Rem                for &Posit<$TYPE> { type Output = Posit<$TYPE> ; fn rem(self, rhs:  Self      ) -> Self::Output { self.quotient_remainder(*rhs).map_or( Posit::<$TYPE>::NAR  , |q_r| q_r.rem) } }
  impl Rem<Posit<$TYPE>>  for &Posit<$TYPE> { type Output = Posit<$TYPE> ; fn rem(self, rhs:  Posit<$TYPE>) -> Self::Output { self.quotient_remainder( rhs).map_or( Posit::<$TYPE>::NAR, |q_r| q_r.rem) } }
  impl Rem<&Posit<$TYPE>> for  Posit<$TYPE> { type Output = Self         ; fn rem(self, rhs: &Posit<$TYPE>) -> Self::Output { self.quotient_remainder(*rhs).map_or( Self::NAR          , |q_r| q_r.rem) } }
  
  impl RemAssign for Posit<$TYPE> { fn rem_assign(&mut self, rhs: Self) { *self = self.rem(rhs) } }
  impl RemAssign<&Self> for Posit<$TYPE> { fn rem_assign(&mut self, rhs: &Self) { *self = self.rem(*rhs) } }
)*};}
impl_operator_traits!{i8 i16 i32 i64 i128 isize}
