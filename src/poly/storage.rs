use memmap2::MmapMut;

use super::Poly;
use std::{fmt::Debug, marker::PhantomData, ops::Range};

// #[derive(Debug, Clone, Copy, PartialEq, Eq)]
// pub enum CompState {
//     Compressed,
//     Decompressed,
// }

pub trait Storage {
    type T: Debug + Clone;
    // type P: Poly;
    fn get(&self, index: Range<usize>) -> &[Self::T];
    fn get_mut(&mut self, index: Range<usize>) -> &mut [Self::T];
    // fn compress_state(&self) -> CompState;
    // fn compress(&mut self);
    // fn decompress(&mut self);
}

// pub struct InMemoryPolyStorage<T: Debug + Clone>(Vec<T>);
// }

impl<T: Debug + Clone> Storage for Vec<T> {
    type T = T;
    fn get(&self, index: Range<usize>) -> &[T] {
        &self[index]
    }

    fn get_mut(&mut self, index: Range<usize>) -> &mut [T] {
        &mut self[index]
    }
}

pub struct MmapPolyStorage<T: Debug + Clone + From<Vec<u8>>> {
    pub mmap: MmapMut,
    pub poly_bytes_len: usize,
    _t: PhantomData<T>,
}

impl<T: Debug + Clone + From<Vec<u8>>> Storage for MmapPolyStorage<T> {
    type T = T;
    fn get(&self, index: Range<usize>) -> &[T] {
        let start = self.poly_bytes_len * index.start;
        let end = self.poly_bytes_len * index.end;
        let bytes = self.mmap.get(start..end).unwrap();
        let polys = bytes
            .chunks_exact(self.poly_bytes_len)
            .map(|bytes| T::from(bytes.to_vec()))
            .collect::<Vec<T>>();
        &polys
    }

    fn get_mut(&mut self, index: Range<usize>) -> &mut [T] {
        todo!()
    }
}
