use crate::{
    poly::{Poly, PolyMatrix},
    utils::{block_size, debug_mem, log_mem},
};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{
    path::{Path, PathBuf},
    sync::Arc,
};
use tokio::{
    runtime::Handle,
    sync::{mpsc, oneshot},
};

#[derive(Debug)]
pub struct SerializedBlock {
    filename: String,
    data: Vec<u8>,
}

#[derive(Debug)]
pub struct SerializedMatrix {
    pub id: String,
    pub blocks: Vec<SerializedBlock>,
}

/// Storage task sent to the I/O worker thread
struct StorageTask {
    serialized_data: Vec<u8>,
    path: PathBuf,
    id: String,
    completion_sender: oneshot::Sender<Result<(), std::io::Error>>,
}

/// Global storage service for handling all matrix I/O operations
pub struct StorageService {
    task_sender: mpsc::UnboundedSender<StorageTask>,
    _worker_handle: tokio::task::JoinHandle<()>,
}

impl StorageService {
    /// Creates a new storage service with configurable concurrent writes
    pub fn new(max_concurrent_writes: usize) -> Self {
        let (task_sender, mut task_receiver) = mpsc::unbounded_channel::<StorageTask>();

        let _worker_handle = tokio::spawn(async move {
            let semaphore = Arc::new(tokio::sync::Semaphore::new(max_concurrent_writes));

            while let Some(task) = task_receiver.recv().await {
                let semaphore = semaphore.clone();

                // Spawn individual I/O tasks up to the concurrency limit
                tokio::spawn(async move {
                    let _permit =
                        semaphore.acquire().await.expect("Semaphore should not be closed");

                    debug_mem(format!("Writing {} to disk", task.id));
                    let result = tokio::fs::write(&task.path, &task.serialized_data).await;

                    match result {
                        Ok(()) => {
                            log_mem(format!("Successfully stored {}", task.id));
                            let _ = task.completion_sender.send(Ok(()));
                        }
                        Err(e) => {
                            eprintln!("Failed to write {}: {}", task.id, e);
                            let _ = task.completion_sender.send(Err(e));
                        }
                    }
                });
            }
        });

        Self { task_sender, _worker_handle }
    }

    /// Queue data to be written to disk - returns immediately
    pub async fn store_data(
        &self,
        data: Vec<u8>,
        path: PathBuf,
        id: String,
    ) -> oneshot::Receiver<Result<(), std::io::Error>> {
        let (completion_sender, completion_receiver) = oneshot::channel();

        let task = StorageTask { serialized_data: data, path, id, completion_sender };

        // Send to I/O worker - this never blocks
        if self.task_sender.send(task).is_err() {
            eprintln!("Storage service is shut down");
        }

        completion_receiver
    }
}

// Global storage service instance using std::sync::OnceLock for thread-safe initialization
static STORAGE_SERVICE: std::sync::OnceLock<StorageService> = std::sync::OnceLock::new();

/// Get the global storage service instance
fn get_storage_service() -> &'static StorageService {
    STORAGE_SERVICE.get_or_init(|| StorageService::new(8)) // Default to 8 concurrent writes
}

/// Handle for managing CPU preprocessing and I/O completion
pub struct StorageHandle {
    cpu_handle: tokio::task::JoinHandle<()>,
    io_completion: oneshot::Receiver<Result<(), std::io::Error>>,
}

impl StorageHandle {
    /// Wait only for CPU preprocessing to complete, return I/O completion handle
    pub async fn wait_cpu_complete(
        self,
    ) -> Result<oneshot::Receiver<Result<(), std::io::Error>>, tokio::task::JoinError> {
        self.cpu_handle.await?;
        Ok(self.io_completion)
    }

    /// Wait for both CPU preprocessing and I/O to complete
    pub async fn wait_all_complete(self) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        self.cpu_handle.await?;
        self.io_completion.await??;
        Ok(())
    }
}

/// CPU-intensive preprocessing function that serializes matrix blocks in parallel
fn preprocess_matrix_for_storage<M>(matrix: M, id: &str) -> SerializedMatrix
where
    M: PolyMatrix + Send + 'static,
{
    let block_size_val = block_size();
    debug_mem(format!("CPU preprocessing started: {id}"));

    // CPU-HEAVY: Extract matrix blocks and serialize in parallel
    let (nrow, ncol) = matrix.size();
    let mut serialized_blocks = Vec::new();

    // Calculate block ranges similar to the original write_to_files implementation
    #[cfg(feature = "disk")]
    let (row_offsets, col_offsets) = {
        use crate::poly::dcrt::matrix::base::disk::block_offsets;
        block_offsets(0..nrow, 0..ncol)
    };
    #[cfg(not(feature = "disk"))]
    let (row_offsets, col_offsets) = (vec![0, nrow], vec![0, ncol]);

    // Process each block in parallel
    use itertools::Itertools;
    let row_windows: Vec<_> = row_offsets.into_iter().tuple_windows().collect();

    for (cur_block_row_idx, next_block_row_idx) in row_windows {
        let col_windows: Vec<_> = col_offsets.clone().into_iter().tuple_windows().collect();

        for (cur_block_col_idx, next_block_col_idx) in col_windows {
            let row_range = cur_block_row_idx..next_block_row_idx;
            let col_range = cur_block_col_idx..next_block_col_idx;

            // CPU-heavy: extract block entries using parallel processing
            let entries = matrix.block_entries(row_range.clone(), col_range.clone());

            // CPU-heavy: serialize each polynomial in the block to compact bytes
            let entries_bytes: Vec<Vec<Vec<u8>>> = entries
                .par_iter()
                .map(|row| row.par_iter().map(|poly| poly.to_compact_bytes()).collect())
                .collect();

            // CPU-heavy: bincode serialization
            let data = bincode::encode_to_vec(&entries_bytes, bincode::config::standard())
                .expect("Failed to serialize matrix block");

            let filename = format!(
                "{}_{}_{}.{}_{}.{}.matrix",
                id, block_size_val, row_range.start, row_range.end, col_range.start, col_range.end
            );

            serialized_blocks.push(SerializedBlock { filename, data });
        }
    }

    debug_mem(format!("CPU preprocessing completed: {} ({} blocks)", id, serialized_blocks.len()));
    SerializedMatrix { id: id.to_string(), blocks: serialized_blocks }
}

/// immediate CPU preprocessing and background I/O
pub fn store_and_drop_matrix<M>(matrix: M, dir: &Path, id: &str) -> StorageHandle
where
    M: PolyMatrix + Send + 'static,
{
    let dir = dir.to_path_buf();
    let id = id.to_owned();

    // IMMEDIATE CPU preprocessing - starts RIGHT NOW
    let (io_completion_sender, io_completion_receiver) = oneshot::channel();

    let cpu_handle = tokio::task::spawn_blocking(move || {
        // CPU-HEAVY: Serialize matrix blocks in parallel
        let serialized_matrix = preprocess_matrix_for_storage(matrix, &id);

        // Early drop after serialization
        debug_mem(format!("Matrix {id} dropped after preprocessing"));

        // Send ALL blocks to I/O service and wait for completion
        let storage_service = get_storage_service();
        let mut io_tasks = Vec::new();

        for block in serialized_matrix.blocks {
            let path = dir.join(&block.filename);
            let block_id = format!("{}::{}", id, block.filename);
            let completion_receiver =
                Handle::current().block_on(storage_service.store_data(block.data, path, block_id));
            io_tasks.push(completion_receiver);
        }

        // Wait for all I/O operations for this matrix to complete
        let final_result = Handle::current().block_on(async {
            let mut all_success = true;
            for io_task in io_tasks {
                match io_task.await {
                    Ok(Ok(())) => {}
                    Ok(Err(e)) => {
                        eprintln!("I/O failed for {id}: {e:?}");
                        all_success = false;
                    }
                    Err(_) => {
                        eprintln!("I/O task cancelled for {id}");
                        all_success = false;
                    }
                }
            }
            if all_success {
                Ok(())
            } else {
                Err(std::io::Error::other("Some I/O operations failed"))
            }
        });

        let _ = io_completion_sender.send(final_result);
    });

    StorageHandle { cpu_handle, io_completion: io_completion_receiver }
}

/// Convert StorageHandle to JoinHandle<()> for backward compatibility
pub fn storage_handle_to_join_handle(storage_handle: StorageHandle) -> tokio::task::JoinHandle<()> {
    tokio::spawn(async move {
        if let Err(e) = storage_handle.wait_all_complete().await {
            eprintln!("Storage operation failed: {e:?}");
        }
    })
}

// // Legacy function for backward compatibility - will be deprecated
// pub fn store_and_drop_matrix_legacy<M>(
//     matrix: M,
//     dir: &Path,
//     id: &str,
// ) -> tokio::task::JoinHandle<()>
// where
//     M: PolyMatrix + Send + 'static,
// {
//     let dir = dir.to_path_buf();
//     let id = id.to_owned();

//     tokio::task::spawn_blocking(move || {
//         log_mem(format!("Storing {id} (legacy)"));
//         Handle::current().block_on(async {
//             matrix.write_to_files(&dir, &id).await;
//         });
//         drop(matrix);
//         log_mem(format!("Stored {id} (legacy)"));
//     })
// }

#[cfg(feature = "debug")]
pub fn store_and_drop_poly<P: Poly + 'static>(
    poly: P,
    dir: &Path,
    id: &str,
) -> tokio::task::JoinHandle<()> {
    let dir = dir.to_path_buf();
    let id = id.to_owned();

    tokio::spawn(async move {
        log_mem(format!("Storing {id}"));
        poly.write_to_file(&dir, &id).await;
        drop(poly);
        log_mem(format!("Stored {id}"));
    })
}
