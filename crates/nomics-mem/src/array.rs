//! Typed numeric arrays backed by HostBuf, with the same allocation policy.
use std::io;
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};

use crate::HostBuf;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for u64 {}
    impl Sealed for u32 {}
    impl Sealed for usize {}
}

/// Numeric types for which every bit pattern, including all zeroes, is valid.
/// Sealed so callers cannot add types with stronger validity/alignment rules.
pub trait HostElement: sealed::Sealed + Copy + std::fmt::Debug {}
impl HostElement for f64 {}
impl HostElement for u64 {}
impl HostElement for u32 {}
impl HostElement for usize {}

/// Fixed-length numeric storage. Allocates through HostBuf before first touch;
/// no temporary Vec is needed when collecting an exact-size iterator.
pub struct HostArray<T: HostElement> {
    buffer: HostBuf,
    element: PhantomData<T>,
}

impl<T: HostElement> HostArray<T> {
    pub fn zeroed(len: usize) -> io::Result<Self> {
        let bytes = len
            .checked_mul(std::mem::size_of::<T>())
            .filter(|&bytes| bytes <= isize::MAX as usize)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "numeric array size overflow")
            })?;
        Ok(Self {
            buffer: HostBuf::new(bytes)?,
            element: PhantomData,
        })
    }

    pub fn collect_exact(values: impl ExactSizeIterator<Item = T>) -> io::Result<Self> {
        let mut array = Self::zeroed(values.len())?;
        for (slot, value) in array.iter_mut().zip(values) {
            *slot = value;
        }
        Ok(array)
    }
}

impl<T: HostElement> Deref for HostArray<T> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        // SAFETY: HostBuf is 2 MiB aligned, the length was checked in zeroed,
        // and the sealed element types admit all bit patterns and have no drop.
        unsafe {
            std::slice::from_raw_parts(
                self.buffer.as_ptr().cast(),
                self.buffer.len() / std::mem::size_of::<T>(),
            )
        }
    }
}

impl<T: HostElement> DerefMut for HostArray<T> {
    fn deref_mut(&mut self) -> &mut [T] {
        // SAFETY: same layout as Deref, with exclusive access to the buffer.
        unsafe {
            std::slice::from_raw_parts_mut(
                self.buffer.as_mut_ptr().cast(),
                self.buffer.len() / std::mem::size_of::<T>(),
            )
        }
    }
}

impl<T: HostElement> std::fmt::Debug for HostArray<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.deref().fmt(f)
    }
}

impl<'a, T: HostElement> IntoIterator for &'a HostArray<T> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a, T: HostElement> IntoIterator for &'a mut HostArray<T> {
    type Item = &'a mut T;
    type IntoIter = std::slice::IterMut<'a, T>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numeric_storage_is_aligned_zeroed_and_mutable() {
        let mut a = HostArray::<f64>::zeroed(300_000).unwrap();
        assert_eq!(a.as_ptr() as usize % crate::HUGE_PAGE, 0);
        assert!(a.iter().all(|&x| x == 0.0));
        a[299_999] = -1.25;
        assert_eq!(a[299_999].to_bits(), (-1.25f64).to_bits());
        let b = HostArray::collect_exact([3u64, 7, 11].into_iter()).unwrap();
        assert_eq!(&*b, &[3, 7, 11]);
        assert!(HostArray::<u32>::zeroed(0).unwrap().is_empty());
        assert!(HostArray::<f64>::zeroed(usize::MAX).is_err());
    }
}
