use std::path::{Path, PathBuf};

pub struct FileHandler {
    folder_path: PathBuf,
}

impl FileHandler {
    pub fn new(folder_path: impl Into<PathBuf>) -> Self {
        Self { folder_path: folder_path.into() }
    }

    pub fn store(&self, data: &[u8], index: usize) {}

    pub fn load(&self, data: &[u8], index: usize) {}
}
