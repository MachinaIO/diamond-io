use crate::{
    poly::{Poly, PolyMatrix},
    utils::{block_size, debug_mem, log_mem},
};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{
    path::{Path, PathBuf},
    time::Instant,
};

// #[derive(Debug)]
// pub struct SerializedBlock {
//     filename: String,
//     data: Vec<u8>,
// }

// #[derive(Debug)]
// pub struct SerializedMatrix {
//     pub id: String,
//     pub blocks: Vec<SerializedBlock>,
// }

// /// Storage task sent to the I/O worker thread
// struct StorageTask {
//     serialized_data: Vec<u8>,
//     path: PathBuf,
//     id: String,
//     completion_sender: oneshot::Sender<Result<(), std::io::Error>>,
// }

// /// Global storage service for handling all matrix I/O operations
// pub struct StorageService {
//     task_sender: mpsc::UnboundedSender<StorageTask>,
//     _worker_handle: tokio::task::JoinHandle<()>,
// }

// impl StorageService {
//     /// Creates a new storage service with configurable concurrent writes
//     pub fn new(max_concurrent_writes: usize) -> Self {
//         let (task_sender, mut task_receiver) = mpsc::unbounded_channel::<StorageTask>();

//         let _worker_handle = tokio::spawn(async move {
//             let semaphore = Arc::new(tokio::sync::Semaphore::new(max_concurrent_writes));

//             while let Some(task) = task_receiver.recv().await {
//                 let semaphore = semaphore.clone();

//                 // Spawn individual I/O tasks up to the concurrency limit
//                 tokio::spawn(async move {
//                     let _permit =
//                         semaphore.acquire().await.expect("Semaphore should not be closed");

//                     debug_mem(format!("Writing {} to disk", task.id));
//                     let result = tokio::fs::write(&task.path, &task.serialized_data).await;

//                     match result {
//                         Ok(()) => {
//                             log_mem(format!("Successfully stored {}", task.id));
//                             let _ = task.completion_sender.send(Ok(()));
//                         }
//                         Err(e) => {
//                             eprintln!("Failed to write {}: {}", task.id, e);
//                             let _ = task.completion_sender.send(Err(e));
//                         }
//                     }
//                 });
//             }
//         });

//         Self { task_sender, _worker_handle }
//     }

//     /// Queue data to be written to disk - returns immediately
//     pub fn store_data(
//         &self,
//         data: Vec<u8>,
//         path: PathBuf,
//         id: String,
//     ) -> oneshot::Receiver<Result<(), std::io::Error>> {
//         let (completion_sender, completion_receiver) = oneshot::channel();

//         let task = StorageTask { serialized_data: data, path, id, completion_sender };

//         // Send to I/O worker - this never blocks
//         if self.task_sender.send(task).is_err() {
//             eprintln!("Storage service is shut down");
//         }

//         completion_receiver
//     }
// }

// // Global storage service instance using std::sync::OnceLock for thread-safe initialization
// static STORAGE_SERVICE: std::sync::OnceLock<StorageService> = std::sync::OnceLock::new();

// /// Get the global storage service instance
// fn get_storage_service() -> &'static StorageService {
//     STORAGE_SERVICE.get_or_init(|| StorageService::new(8)) // Default to 8 concurrent writes
// }

// /// Handle for managing CPU preprocessing and I/O completion
// pub struct StorageHandle {
//     cpu_handle: tokio::task::JoinHandle<()>,
//     io_completion: oneshot::Receiver<Result<(), std::io::Error>>,
// }

// impl StorageHandle {
//     /// Wait only for CPU preprocessing to complete, return I/O completion handle
//     pub async fn wait_cpu_complete(
//         self,
//     ) -> Result<oneshot::Receiver<Result<(), std::io::Error>>, tokio::task::JoinError> {
//         self.cpu_handle.await?;
//         Ok(self.io_completion)
//     }

//     /// Wait for both CPU preprocessing and I/O to complete
//     pub async fn wait_all_complete(self) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
// {         self.cpu_handle.await?;
//         self.io_completion.await??;
//         Ok(())
//     }
// }

// /// CPU-intensive preprocessing function that serializes matrix blocks in parallel
// fn preprocess_matrix_for_storage<M>(matrix: M, id: &str) -> SerializedMatrix
// where
//     M: PolyMatrix + Send + 'static,
// {
//     let start = Instant::now();
//     let block_size_val = block_size();
//     debug_mem(format!("CPU preprocessing started: {id}"));

//     // CPU-HEAVY: Extract matrix blocks and serialize in parallel
//     let (nrow, ncol) = matrix.size();
//     let mut serialized_blocks = Vec::new();

//     // Calculate block ranges similar to the original write_to_files implementation
//     #[cfg(feature = "disk")]
//     let (row_offsets, col_offsets) = {
//         use crate::poly::dcrt::matrix::base::disk::block_offsets;
//         block_offsets(0..nrow, 0..ncol)
//     };
//     #[cfg(not(feature = "disk"))]
//     let (row_offsets, col_offsets) = (vec![0, nrow], vec![0, ncol]);

//     // Process each block in parallel
//     use itertools::Itertools;
//     let row_windows: Vec<_> = row_offsets.into_iter().tuple_windows().collect();

//     for (cur_block_row_idx, next_block_row_idx) in row_windows {
//         let col_windows: Vec<_> = col_offsets.clone().into_iter().tuple_windows().collect();

//         for (cur_block_col_idx, next_block_col_idx) in col_windows {
//             let row_range = cur_block_row_idx..next_block_row_idx;
//             let col_range = cur_block_col_idx..next_block_col_idx;

//             // CPU-heavy: extract block entries using parallel processing
//             let entries = matrix.block_entries(row_range.clone(), col_range.clone());

//             // CPU-heavy: serialize each polynomial in the block to compact bytes
//             let entries_bytes: Vec<Vec<Vec<u8>>> = entries
//                 .par_iter()
//                 .map(|row| row.par_iter().map(|poly| poly.to_compact_bytes()).collect())
//                 .collect();

//             // CPU-heavy: bincode serialization
//             let data = bincode::encode_to_vec(&entries_bytes, bincode::config::standard())
//                 .expect("Failed to serialize matrix block");

//             let filename = format!(
//                 "{}_{}_{}.{}_{}.{}.matrix",
//                 id, block_size_val, row_range.start, row_range.end, col_range.start,
// col_range.end             );
//             info!("SerializedBlock len {}", data.len());

//             serialized_blocks.push(SerializedBlock { filename, data });
//         }
//     }
//     let elapsed = start.elapsed();
//     log_mem(format!(
//         "CPU preprocessing completed: {} ({} blocks) {elapsed:?}",
//         id,
//         serialized_blocks.len()
//     ));
//     SerializedMatrix { id: id.to_string(), blocks: serialized_blocks }
// }

// /// immediate CPU preprocessing and background I/O
// pub fn store_and_drop_matrix<M>(matrix: M, dir: &Path, id: &str) -> StorageHandle
// where
//     M: PolyMatrix + Send + 'static,
// {
//     let dir = dir.to_path_buf();
//     let id = id.to_owned();

//     // IMMEDIATE CPU preprocessing - starts RIGHT NOW
//     let (io_completion_sender, io_completion_receiver) = oneshot::channel();

//     let cpu_handle = tokio::task::spawn_blocking(move || {
//         let serialized_matrix = preprocess_matrix_for_storage(matrix, &id);
//         debug_mem(format!("Matrix {id} dropped after preprocessing"));
//         let storage_service = get_storage_service();
//         let mut io_tasks = Vec::new();

//         for block in serialized_matrix.blocks {
//             let path = dir.join(&block.filename);
//             let block_id = format!("{}::{}", id, block.filename);
//             let completion_receiver = storage_service.store_data(block.data, path, block_id);
//             io_tasks.push(completion_receiver);
//         }

//         // Wait for all I/O operations for this matrix to complete
//         let final_result = Handle::current().block_on(async {
//             let mut all_success = true;
//             for io_task in io_tasks {
//                 match io_task.await {
//                     Ok(Ok(())) => {}
//                     Ok(Err(e)) => {
//                         eprintln!("I/O failed for {id}: {e:?}");
//                         all_success = false;
//                     }
//                     Err(_) => {
//                         eprintln!("I/O task cancelled for {id}");
//                         all_success = false;
//                     }
//                 }
//             }
//             if all_success {
//                 Ok(())
//             } else {
//                 Err(std::io::Error::other("Some I/O operations failed"))
//             }
//         });

//         let _ = io_completion_sender.send(final_result);
//     });

//     StorageHandle { cpu_handle, io_completion: io_completion_receiver }
// }

// /// Convert StorageHandle to JoinHandle<()> for backward compatibility
// pub fn storage_handle_to_join_handle(storage_handle: StorageHandle) ->
// tokio::task::JoinHandle<()> {     tokio::spawn(async move {
//         if let Err(e) = storage_handle.wait_all_complete().await {
//             eprintln!("Storage operation failed: {e:?}");
//         }
//     })
// }

/// Optimized streaming storage - processes and writes blocks incrementally
/// This version reduces peak memory usage by streaming blocks to disk
pub fn store_and_drop_matrix_streaming<M>(
    matrix: M,
    dir: &Path,
    id: &str,
) -> tokio::task::JoinHandle<()>
where
    M: PolyMatrix + Send + 'static,
{
    let dir = dir.to_path_buf();
    let id = id.to_owned();

    tokio::task::spawn_blocking(move || {
        let start = Instant::now();
        debug_mem(format!("Streaming storage started: {id}"));

        // Use larger buffer to prevent serialization blocking on I/O
        let buffer_size =
            std::env::var("STORAGE_BUFFER_SIZE").ok().and_then(|s| s.parse().ok()).unwrap_or(20);
        let (sender, receiver) = std::sync::mpsc::sync_channel::<(PathBuf, Vec<u8>)>(buffer_size);

        // Track buffer memory usage
        let buffer_memory = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let buffer_memory_io = buffer_memory.clone();

        // Simple I/O thread
        let io_handle = std::thread::spawn(move || {
            let mut blocks_written = 0;
            let mut total_bytes_written = 0;
            for (path, data) in receiver {
                let data_len = data.len();
                if let Err(e) = std::fs::write(&path, &data) {
                    eprintln!("Failed to write {}: {}", path.display(), e);
                } else {
                    blocks_written += 1;
                    total_bytes_written += data_len;
                    // Subtract from buffer memory after write
                    buffer_memory_io.fetch_sub(data_len, std::sync::atomic::Ordering::Relaxed);
                    log_mem(format!(
                        "Block {} written ({} bytes, {} MB total, buffer: {} MB)",
                        blocks_written,
                        data_len,
                        total_bytes_written / 1_000_000,
                        buffer_memory_io.load(std::sync::atomic::Ordering::Relaxed) / 1_000_000
                    ));
                }
            }
            log_mem(format!(
                "I/O thread completed: {} blocks written, {} MB total",
                blocks_written,
                total_bytes_written / 1_000_000
            ));
        });

        // Process and write block-by-block
        stream_matrix_blocks(&matrix, &dir, &id, sender, buffer_memory);

        io_handle.join().unwrap();
        drop(matrix);

        let elapsed = start.elapsed();
        log_mem(format!("Streaming storage completed: {id} in {elapsed:?}"));
    })
}

/// Stream matrix blocks for immediate processing and I/O
fn stream_matrix_blocks<M>(
    matrix: &M,
    dir: &Path,
    id: &str,
    sender: std::sync::mpsc::SyncSender<(PathBuf, Vec<u8>)>,
    buffer_memory: std::sync::Arc<std::sync::atomic::AtomicUsize>,
) where
    M: PolyMatrix,
{
    let block_size_val = block_size();
    let (nrow, ncol) = matrix.size();
    let block_start = Instant::now();

    #[cfg(feature = "disk")]
    let (row_offsets, col_offsets) = {
        use crate::poly::dcrt::matrix::base::disk::block_offsets;
        block_offsets(0..nrow, 0..ncol)
    };
    #[cfg(not(feature = "disk"))]
    let (row_offsets, col_offsets) = (vec![0, nrow], vec![0, ncol]);

    let mut blocks_processed = 0;
    let mut total_serialized_bytes = 0;

    // STREAMING: Process each block immediately
    use itertools::Itertools;
    let row_windows: Vec<_> = row_offsets.into_iter().tuple_windows().collect();
    let total_blocks = row_windows.len();

    for (cur_block_row_idx, next_block_row_idx) in &row_windows {
        let col_windows: Vec<_> = col_offsets.clone().into_iter().tuple_windows().collect();

        for (cur_block_col_idx, next_block_col_idx) in &col_windows {
            let row_range = *cur_block_row_idx..*next_block_row_idx;
            let col_range = *cur_block_col_idx..*next_block_col_idx;

            let block_cpu_start = Instant::now();

            // CPU work for this block only
            let entries = matrix.block_entries(row_range.clone(), col_range.clone());
            let entries_bytes: Vec<Vec<Vec<u8>>> = entries
                .par_iter()
                .map(|row| row.par_iter().map(|poly| poly.to_compact_bytes()).collect())
                .collect();

            let data = bincode::encode_to_vec(&entries_bytes, bincode::config::standard())
                .expect("Failed to serialize matrix block");

            let filename = format!(
                "{}_{}_{}.{}_{}.{}.matrix",
                id, block_size_val, row_range.start, row_range.end, col_range.start, col_range.end
            );

            let block_cpu_elapsed = block_cpu_start.elapsed();
            blocks_processed += 1;
            total_serialized_bytes += data.len();

            let data_len = data.len();
            debug_mem(format!(
                "Block {}/{} processed in {block_cpu_elapsed:?} ({} bytes)",
                blocks_processed,
                total_blocks * col_windows.len(),
                data_len
            ));

            buffer_memory.fetch_add(data_len, std::sync::atomic::Ordering::Relaxed);
            let current_buffer = buffer_memory.load(std::sync::atomic::Ordering::Relaxed);
            if current_buffer > 10_000_000_000 {
                // Warn if buffer > 10GB
                log_mem(format!(
                    "WARNING: Buffer memory high: {} GB (block {})",
                    current_buffer / 1_000_000_000,
                    blocks_processed
                ));
            }

            // Send to I/O
            if sender.send((dir.join(filename), data)).is_err() {
                eprintln!("I/O thread died, stopping processing");
                break;
            }
        }
    }

    let total_elapsed = block_start.elapsed();
    log_mem(format!(
        "Matrix streaming completed: {id} - {} blocks, {} total bytes in {total_elapsed:?}",
        blocks_processed, total_serialized_bytes
    ));
}

/// Read a matrix that was split into blocks back into a single matrix.
pub fn read_matrix_from_blocks<M>(
    params: &<M::P as Poly>::Params,
    nrow: usize,
    ncol: usize,
    dir: &Path,
    id: &str,
) -> M
where
    M: PolyMatrix,
{
    let start = Instant::now();
    debug_mem(format!("Reading matrix from blocks: {id}"));
    let matrix = M::read_from_files(params, nrow, ncol, dir, id);

    let elapsed = start.elapsed();
    log_mem(format!("Matrix read completed: {id} in {elapsed:?}"));

    matrix
}

/// Read matrix metadata to determine dimensions from block files.
pub fn read_matrix_metadata(dir: &Path, id: &str) -> Result<(usize, usize), std::io::Error> {
    use std::fs;

    let block_size_val = block_size();

    // Find all files matching the pattern
    let pattern = format!("{}_{}_{}", id, block_size_val, "");
    let mut max_row = 0;
    let mut max_col = 0;

    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let filename = entry.file_name();
        let filename_str = filename.to_string_lossy();

        if filename_str.starts_with(&pattern) && filename_str.ends_with(".matrix") {
            // Parse the filename to extract ranges
            // Format: id_blocksize_rowstart.rowend_colstart.colend.matrix
            // First split by underscore to get [id, blocksize, ranges...]
            let parts: Vec<&str> = filename_str.split('_').collect();
            if parts.len() >= 3 {
                // The range part is everything after id_blocksize_
                let range_part = &filename_str[pattern.len()..];

                // Now we have something like "0.10_0.15.matrix"
                // Split by underscore to get row and column parts
                let range_parts: Vec<&str> =
                    range_part.trim_end_matches(".matrix").split('_').collect();

                if range_parts.len() >= 2 {
                    // Extract row end from "0.10"
                    if let Some(row_end) = range_parts[0].split('.').nth(1) {
                        if let Ok(row) = row_end.parse::<usize>() {
                            max_row = max_row.max(row);
                        }
                    }

                    // Extract col end from "0.15"
                    if let Some(col_end) = range_parts[1].split('.').nth(1) {
                        if let Ok(col) = col_end.parse::<usize>() {
                            max_col = max_col.max(col);
                        }
                    }
                }
            }
        }
    }

    if max_row == 0 || max_col == 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("No matrix blocks found for id: {}", id),
        ));
    }

    Ok((max_row, max_col))
}

pub fn read_matrix_streaming<M, F>(
    params: &<M::P as Poly>::Params,
    nrow: usize,
    ncol: usize,
    dir: &Path,
    id: &str,
    mut process_block: F,
) -> Result<(), std::io::Error>
where
    M: PolyMatrix,
    F: FnMut(usize, usize, Vec<Vec<M::P>>),
{
    let block_size_val = block_size();

    #[cfg(feature = "disk")]
    let (row_offsets, col_offsets) = {
        use crate::poly::dcrt::matrix::base::disk::block_offsets;
        block_offsets(0..nrow, 0..ncol)
    };
    #[cfg(not(feature = "disk"))]
    let (row_offsets, col_offsets) = (vec![0, nrow], vec![0, ncol]);

    use itertools::Itertools;
    let row_windows: Vec<_> = row_offsets.into_iter().tuple_windows().collect();

    for (row_idx, (cur_block_row_idx, next_block_row_idx)) in row_windows.iter().enumerate() {
        let col_windows: Vec<_> = col_offsets.clone().into_iter().tuple_windows().collect();

        for (col_idx, (cur_block_col_idx, next_block_col_idx)) in col_windows.iter().enumerate() {
            let row_range = *cur_block_row_idx..*next_block_row_idx;
            let col_range = *cur_block_col_idx..*next_block_col_idx;

            let filename = format!(
                "{}_{}_{}.{}_{}.{}.matrix",
                id, block_size_val, row_range.start, row_range.end, col_range.start, col_range.end
            );

            let path = dir.join(&filename);
            let bytes = std::fs::read(&path)?;

            let entries_bytes: Vec<Vec<Vec<u8>>> =
                bincode::decode_from_slice(&bytes, bincode::config::standard())
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?
                    .0;

            // Deserialize the block
            let block_entries: Vec<Vec<M::P>> = entries_bytes
                .par_iter()
                .map(|row| {
                    row.par_iter()
                        .map(|entry_bytes| M::P::from_compact_bytes(params, entry_bytes))
                        .collect()
                })
                .collect();

            // Process this block
            process_block(row_idx, col_idx, block_entries);
        }
    }

    Ok(())
}

#[cfg(feature = "debug")]
pub fn store_and_drop_poly<P: Poly>(poly: P, dir: &Path, id: &str) {
    log_mem(format!("Storing {id}"));
    poly.write_to_file(dir, id);
    drop(poly);
    log_mem(format!("Stored {id}"));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{
        dcrt::{DCRTPoly, DCRTPolyMatrix, DCRTPolyParams},
        Poly,
    };
    use tempfile::tempdir;

    #[tokio::test]
    async fn test_matrix_write_read_roundtrip() {
        // Create test parameters
        let params = DCRTPolyParams::new(4, 2, 17, 1);

        // Create a small test matrix
        let nrow = 4;
        let ncol = 4;
        let mut matrix = DCRTPolyMatrix::zero(&params, nrow, ncol);

        // Fill with test data
        for i in 0..nrow {
            for j in 0..ncol {
                let poly = DCRTPoly::const_int(&params, i * ncol + j);
                matrix.set_entry(i, j, poly);
            }
        }

        // Create temp directory
        let temp_dir = tempdir().unwrap();
        let dir_path = temp_dir.path();
        let matrix_id = "test_matrix";

        // Write matrix using streaming storage
        let write_handle = store_and_drop_matrix_streaming(matrix.clone(), dir_path, matrix_id);
        write_handle.await.unwrap();

        // Read matrix back
        let read_matrix =
            read_matrix_from_blocks::<DCRTPolyMatrix>(&params, nrow, ncol, dir_path, matrix_id);

        assert_eq!(matrix, read_matrix);
    }

    #[tokio::test]
    async fn test_matrix_metadata_read() {
        // Create test parameters
        let params = DCRTPolyParams::new(4, 2, 17, 1);

        // Create test matrix with known dimensions
        let nrow = 10;
        let ncol = 15;
        let matrix = DCRTPolyMatrix::zero(&params, nrow, ncol);

        // Create temp directory
        let temp_dir = tempdir().unwrap();
        let dir_path = temp_dir.path();
        let matrix_id = "test_metadata";

        // Write matrix
        let write_handle = store_and_drop_matrix_streaming(matrix, dir_path, matrix_id);
        write_handle.await.unwrap();

        // Read metadata
        let (read_nrow, read_ncol) = read_matrix_metadata(dir_path, matrix_id).unwrap();

        // Verify dimensions
        assert_eq!(read_nrow, nrow);
        assert_eq!(read_ncol, ncol);
    }

    #[tokio::test]
    async fn test_matrix_streaming_read() {
        // Create test parameters
        let params = DCRTPolyParams::new(4, 2, 17, 1);

        // Create larger test matrix to ensure multiple blocks
        let nrow = 168;
        let ncol = 168;
        let mut matrix = DCRTPolyMatrix::zero(&params, nrow, ncol);

        // Fill with test data
        for i in 0..nrow {
            for j in 0..ncol {
                let poly = DCRTPoly::const_int(&params, i * ncol + j);
                matrix.set_entry(i, j, poly);
            }
        }

        // Create temp directory
        let temp_dir = tempdir().unwrap();
        let dir_path = temp_dir.path();
        let matrix_id = "test_streaming";

        // Write matrix
        let write_handle = store_and_drop_matrix_streaming(matrix.clone(), dir_path, matrix_id);
        write_handle.await.unwrap();

        // Read matrix using streaming
        let mut block_count = 0;
        let mut total_elements = 0;

        read_matrix_streaming::<DCRTPolyMatrix, _>(
            &params,
            nrow,
            ncol,
            dir_path,
            matrix_id,
            |_row_idx, _col_idx, block| {
                block_count += 1;
                for row in &block {
                    total_elements += row.len();
                }
            },
        )
        .unwrap();

        // Verify we read all blocks and elements
        assert!(block_count > 0);
        assert_eq!(total_elements, nrow * ncol);
    }
}
