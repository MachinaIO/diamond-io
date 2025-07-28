use crate::{
    poly::{Poly, PolyMatrix},
    utils::{block_size, debug_mem, log_mem},
};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{
    path::Path,
    sync::{Arc, Mutex, OnceLock},
    time::Instant,
};
use tokio::task::JoinHandle;
use tracing::info;

#[derive(Debug)]
pub struct SerializedMatrix {
    pub id: String,
    pub filename: String,
    pub data: Vec<u8>,
}

static WRITE_HANDLES: OnceLock<Arc<Mutex<Vec<JoinHandle<()>>>>> = OnceLock::new();

/// Initialize the storage system
pub fn init_storage_system() {
    WRITE_HANDLES
        .set(Arc::new(Mutex::new(Vec::new())))
        .expect("Storage system already initialized");
}

/// Wait for all pending writes to complete
pub async fn wait_for_all_writes() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let handles = WRITE_HANDLES.get().ok_or("Storage system not initialized")?;

    let handles_vec: Vec<JoinHandle<()>> = {
        let mut guard = handles.lock().unwrap();
        std::mem::take(guard.as_mut())
    };

    log_mem(format!("Waiting for {} pending writes to complete", handles_vec.len()));

    for handle in handles_vec {
        if let Err(e) = handle.await {
            eprintln!("Write task failed: {}", e);
        }
    }

    log_mem("All writes completed");
    Ok(())
}

fn preprocess_matrix_for_storage<M>(matrix: M, id: &str) -> SerializedMatrix
where
    M: PolyMatrix + Send + 'static,
{
    let start = Instant::now();
    let block_size_val = block_size();
    debug_mem(format!("CPU preprocessing started: {id}"));

    let (nrow, ncol) = matrix.size();
    let row_range = 0..nrow;
    let col_range = 0..ncol;

    let entries = matrix.block_entries(row_range.clone(), col_range.clone());
    let entries_bytes: Vec<Vec<Vec<u8>>> = entries
        .par_iter()
        .map(|row| row.par_iter().map(|poly| poly.to_compact_bytes()).collect())
        .collect();
    let data = bincode::encode_to_vec(&entries_bytes, bincode::config::standard())
        .expect("Failed to serialize matrix");

    let filename = format!(
        "{}_{}_{}.{}_{}.{}.matrix",
        id, block_size_val, row_range.start, row_range.end, col_range.start, col_range.end
    );

    let elapsed = start.elapsed();
    info!("Serialized matrix {} len {} bytes in {elapsed:?}", id, data.len());
    log_mem(format!("CPU preprocessing completed: {id} ({} bytes)", data.len()));

    SerializedMatrix { id: id.to_string(), filename, data }
}

/// CPU preprocessing (blocking) + background I/O (non-blocking)
pub fn store_and_drop_matrix<M>(matrix: M, dir: &Path, id: &str)
where
    M: PolyMatrix + Send + 'static,
{
    let dir = dir.to_path_buf();
    let id = id.to_owned();

    let serialized_matrix = preprocess_matrix_for_storage(matrix, &id);
    log_mem(format!("Matrix {id} dropped after preprocessing"));

    // Spawn async I/O task - non-blocking
    let write_handle = tokio::spawn(async move {
        let path = dir.join(&serialized_matrix.filename);

        match tokio::fs::write(&path, &serialized_matrix.data).await {
            Ok(_) => {
                log_mem(format!(
                    "Matrix {} written to {} ({} bytes)",
                    serialized_matrix.id,
                    serialized_matrix.filename,
                    serialized_matrix.data.len()
                ));
            }
            Err(e) => {
                eprintln!("Failed to write {}: {}", path.display(), e);
            }
        }
    });

    // Store handle for later synchronization
    if let Some(handles) = WRITE_HANDLES.get() {
        handles.lock().unwrap().push(write_handle);
    } else {
        eprintln!("Warning: Storage system not initialized, falling back to blocking write");
        // Fallback: block on the async operation
        tokio::task::block_in_place(|| tokio::runtime::Handle::current().block_on(write_handle))
            .unwrap();
    }
}

#[cfg(feature = "debug")]
pub fn store_and_drop_poly<P: Poly>(poly: P, dir: &Path, id: &str) {
    log_mem(format!("Storing {id}"));
    poly.write_to_file(dir, id);
    drop(poly);
    log_mem(format!("Stored {id}"));
}
