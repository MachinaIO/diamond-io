use crate::{
    poly::{Poly, PolyMatrix},
    utils::{block_size, debug_mem, log_mem},
};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicUsize, Ordering},
        mpsc::{sync_channel, Receiver, SyncSender},
        Arc,
    },
    thread,
    time::Instant,
};

/// Optimized configuration for large matrix storage
#[derive(Clone)]
pub struct StorageConfig {
    /// Number of parallel serialization workers
    pub serialization_workers: usize,
    /// Number of I/O writer threads
    pub io_writers: usize,
    /// Buffer size in number of blocks (not fixed count)
    pub buffer_blocks: usize,
    /// Block size override (if None, uses default)
    pub block_size_override: Option<usize>,
    /// Enable direct I/O (bypasses OS cache)
    pub direct_io: bool,
    /// Maximum memory for buffers in bytes
    pub max_buffer_memory: usize,
}

impl Default for StorageConfig {
    fn default() -> Self {
        let num_cpus = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8);
        Self {
            serialization_workers: num_cpus.saturating_sub(2).max(1),
            io_writers: 4,
            buffer_blocks: 100,
            block_size_override: Some(1000),
            direct_io: true,
            max_buffer_memory: 4_000_000_000,
        }
    }
}

/// Represents a serialized block ready for I/O
struct SerializedBlock {
    path: PathBuf,
    data: Vec<u8>,
}

/// Optimized streaming storage with parallel serialization and multiple I/O threads
pub fn store_and_drop_matrix_streaming_optimized<M>(
    matrix: M,
    dir: &Path,
    id: &str,
    config: StorageConfig,
) -> tokio::task::JoinHandle<()>
where
    M: PolyMatrix + Send + Sync + 'static,
{
    let dir = dir.to_path_buf();
    let id = id.to_owned();

    tokio::task::spawn_blocking(move || {
        let start = Instant::now();
        log_mem(format!("Optimized streaming storage started: {id}"));

        // Shared atomic counters for monitoring
        let blocks_serialized = Arc::new(AtomicUsize::new(0));
        let blocks_written = Arc::new(AtomicUsize::new(0));
        let buffer_memory = Arc::new(AtomicUsize::new(0));
        let total_bytes = Arc::new(AtomicUsize::new(0));

        // Create per-thread channels for parallel pipeline
        let buffer_per_writer = config.buffer_blocks / config.io_writers.max(1);
        let (serialize_txs, serialize_rxs): (Vec<_>, Vec<_>) = (0..config.io_writers)
            .map(|_| sync_channel::<SerializedBlock>(buffer_per_writer))
            .unzip();

        // Spawn I/O writer threads with individual channels
        let io_handles = spawn_io_writers_dedicated(
            serialize_rxs,
            Arc::clone(&blocks_written),
            Arc::clone(&buffer_memory),
            Arc::clone(&total_bytes),
            config.direct_io,
        );

        // Perform parallel block serialization with round-robin distribution
        parallel_serialize_matrix_distributed(
            &matrix,
            &dir,
            &id,
            serialize_txs,
            config,
            blocks_serialized,
            buffer_memory,
        );

        // Wait for all I/O to complete
        for handle in io_handles {
            handle.join().expect("I/O thread panicked");
        }

        drop(matrix);
        let elapsed = start.elapsed();
        let total_mb = total_bytes.load(Ordering::Relaxed) / 1_000_000;
        let throughput =
            if elapsed.as_secs() > 0 { total_mb as f64 / elapsed.as_secs_f64() } else { 0.0 };

        log_mem(format!(
            "Optimized streaming completed: {id} - {} MB in {elapsed:?} ({:.1} MB/s)",
            total_mb, throughput
        ));
    })
}

/// Spawn multiple I/O writer threads with dedicated channels
fn spawn_io_writers_dedicated(
    receivers: Vec<Receiver<SerializedBlock>>,
    blocks_written: Arc<AtomicUsize>,
    buffer_memory: Arc<AtomicUsize>,
    total_bytes: Arc<AtomicUsize>,
    direct_io: bool,
) -> Vec<thread::JoinHandle<()>> {
    receivers
        .into_iter()
        .enumerate()
        .map(|(writer_id, receiver)| {
            let blocks_written = Arc::clone(&blocks_written);
            let buffer_memory = Arc::clone(&buffer_memory);
            let total_bytes = Arc::clone(&total_bytes);

            thread::spawn(move || {
                debug_mem(format!("I/O writer {} started", writer_id));

                while let Ok(block) = receiver.recv() {
                    let data_len = block.data.len();

                    // Write with optional direct I/O
                    let write_result = if direct_io {
                        write_direct_io(&block.path, &block.data)
                    } else {
                        std::fs::write(&block.path, &block.data)
                    };

                    match write_result {
                        Ok(()) => {
                            blocks_written.fetch_add(1, Ordering::Relaxed);
                            total_bytes.fetch_add(data_len, Ordering::Relaxed);
                            buffer_memory.fetch_sub(data_len, Ordering::Relaxed);

                            let written = blocks_written.load(Ordering::Relaxed);
                            if written % 10 == 0 {
                                // Log every 10 blocks
                                log_mem(format!(
                                    "Writer {}: {} blocks written, buffer: {} MB",
                                    writer_id,
                                    written,
                                    buffer_memory.load(Ordering::Relaxed) / 1_000_000
                                ));
                            }
                        }
                        Err(e) => {
                            eprintln!(
                                "Writer {} failed to write {}: {}",
                                writer_id,
                                block.path.display(),
                                e
                            );
                        }
                    }
                }

                debug_mem(format!("I/O writer {} completed", writer_id));
            })
        })
        .collect()
}

/// Parallel serialization of matrix blocks with distributed channels
fn parallel_serialize_matrix_distributed<M>(
    matrix: &M,
    dir: &Path,
    id: &str,
    senders: Vec<SyncSender<SerializedBlock>>,
    config: StorageConfig,
    blocks_serialized: Arc<AtomicUsize>,
    buffer_memory: Arc<AtomicUsize>,
) where
    M: PolyMatrix + Sync,
{
    let block_size_val = config.block_size_override.unwrap_or_else(block_size);
    let (nrow, ncol) = matrix.size();

    log_mem(format!(
        "Starting parallel serialization: {} ({}x{}) with block size {}",
        id, nrow, ncol, block_size_val
    ));

    // Calculate block ranges
    #[cfg(feature = "disk")]
    let (row_offsets, col_offsets) = {
        use crate::poly::dcrt::matrix::base::disk::block_offsets_with_size;
        block_offsets_with_size(0..nrow, 0..ncol, block_size_val)
    };
    #[cfg(not(feature = "disk"))]
    let (row_offsets, col_offsets) = {
        use crate::poly::dcrt::matrix::base::memory::block_offsets_with_size;
        block_offsets_with_size(0..nrow, 0..ncol, block_size_val)
    };

    // Prepare all block coordinates
    use itertools::Itertools;
    let block_coords: Vec<_> = row_offsets
        .iter()
        .tuple_windows()
        .flat_map(|(r1, r2)| {
            col_offsets.iter().tuple_windows().map(move |(c1, c2)| (*r1, *r2, *c1, *c2))
        })
        .collect();

    let total_blocks = block_coords.len();
    log_mem(format!("Processing {} blocks in parallel", total_blocks));

    // Round-robin sender selection
    let sender_idx = Arc::new(AtomicUsize::new(0));
    let num_senders = senders.len();

    // Use rayon to process blocks in parallel chunks
    let chunk_size = (total_blocks / config.serialization_workers).max(1);

    block_coords.chunks(chunk_size).for_each(|chunk| {
        chunk.iter().for_each(|&(r1, r2, c1, c2)| {
            let row_range = r1..r2;
            let col_range = c1..c2;

            // Check memory limit before serializing
            while buffer_memory.load(Ordering::Relaxed) > config.max_buffer_memory {
                std::thread::sleep(std::time::Duration::from_millis(10));
            }

            let serialize_start = Instant::now();

            // Extract and serialize block
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

            let data_len = data.len();
            buffer_memory.fetch_add(data_len, Ordering::Relaxed);

            let block = SerializedBlock { path: dir.join(filename), data };

            // Send to I/O threads using round-robin distribution
            let sender_index = sender_idx.fetch_add(1, Ordering::Relaxed) % num_senders;
            if senders[sender_index].send(block).is_err() {
                eprintln!("I/O thread {} died, stopping serialization", sender_index);
                return;
            }

            let serialized = blocks_serialized.fetch_add(1, Ordering::Relaxed) + 1;
            let serialize_elapsed = serialize_start.elapsed();

            if serialized % 10 == 0 || data_len > 500_000_000 {
                log_mem(format!(
                    "Serialized {}/{} blocks ({:.1}%), last block: {} MB in {:?}",
                    serialized,
                    total_blocks,
                    (serialized as f64 / total_blocks as f64) * 100.0,
                    data_len / 1_000_000,
                    serialize_elapsed
                ));
            }
        });
    });
}

/// Write with direct I/O (bypasses OS cache)
fn write_direct_io(path: &Path, data: &[u8]) -> std::io::Result<()> {
    #[cfg(target_os = "linux")]
    {
        use std::{io::Write, os::unix::fs::OpenOptionsExt};

        // Use O_DIRECT flag to bypass page cache
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .custom_flags(libc::O_DIRECT)
            .open(path)?;

        // O_DIRECT requires aligned buffers and sizes
        let aligned_data = align_buffer_for_direct_io(data);
        file.write_all(&aligned_data)?;
        file.sync_all()?; // Ensure data reaches disk
        Ok(())
    }

    #[cfg(target_os = "macos")]
    {
        use std::{
            io::Write,
            os::unix::{fs::OpenOptionsExt, io::AsRawFd},
        };

        // macOS uses F_NOCACHE fcntl instead of O_DIRECT
        let mut file =
            std::fs::OpenOptions::new().write(true).create(true).truncate(true).open(path)?;

        // Set F_NOCACHE to bypass buffer cache
        unsafe {
            let fd = file.as_raw_fd();
            libc::fcntl(fd, libc::F_NOCACHE, 1);
        }

        file.write_all(data)?;
        file.sync_all()?; // Ensure data reaches disk
        Ok(())
    }

    #[cfg(target_os = "windows")]
    {
        use std::{io::Write, os::windows::fs::OpenOptionsExt};

        // Windows uses FILE_FLAG_NO_BUFFERING
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .custom_flags(0x20000000) // FILE_FLAG_NO_BUFFERING
            .open(path)?;

        // Windows direct I/O requires sector-aligned buffers
        let aligned_data = align_buffer_for_direct_io(data);
        file.write_all(&aligned_data)?;
        file.sync_all()?;
        Ok(())
    }

    #[cfg(not(any(target_os = "linux", target_os = "macos", target_os = "windows")))]
    {
        // Fallback for unsupported platforms
        use std::io::Write;
        let mut file = std::fs::File::create(path)?;
        file.write_all(data)?;
        file.sync_all()?;
        Ok(())
    }
}

/// Align buffer for direct I/O requirements
#[cfg(target_os = "linux")]
fn align_buffer_for_direct_io(data: &[u8]) -> Vec<u8> {
    const ALIGNMENT: usize = 4096; // Common page size
    let aligned_size = (data.len() + ALIGNMENT - 1) & !(ALIGNMENT - 1);

    let mut aligned_buffer = vec![0u8; aligned_size];
    aligned_buffer[..data.len()].copy_from_slice(data);
    aligned_buffer
}

/// Public API function that uses optimized storage with default config
pub fn store_matrix_optimized<M>(matrix: M, dir: &Path, id: &str) -> tokio::task::JoinHandle<()>
where
    M: PolyMatrix + Send + Sync + 'static,
{
    store_and_drop_matrix_streaming_optimized(matrix, dir, id, StorageConfig::default())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        poly::dcrt::{DCRTPolyMatrix, DCRTPolyParams},
        utils::{create_random_poly, init_tracing},
    };
    use std::fs;
    use tempfile::TempDir;
    use tokio::runtime::Runtime;

    fn create_test_matrix(rows: usize, cols: usize) -> DCRTPolyMatrix {
        let params = DCRTPolyParams::new(8192, 7, 51, 17);
        let mut matrix = DCRTPolyMatrix::zero(&params, rows, cols);

        for i in 0..rows {
            for j in 0..cols {
                let poly = create_random_poly(&params);
                matrix.set_entry(i, j, poly);
            }
        }
        matrix
    }

    #[test]
    fn test_storage_config_default() {
        let config = StorageConfig::default();
        assert!(config.serialization_workers > 0);
        assert!(config.io_writers > 0);
        assert!(config.buffer_blocks > 0);
        assert!(config.max_buffer_memory > 0);
    }

    #[test]
    fn test_optimized_storage_small_matrix() {
        let rt = Runtime::new().unwrap();
        let temp_dir = TempDir::new().unwrap();
        let matrix = create_test_matrix(10, 10);

        let config = StorageConfig {
            serialization_workers: 2,
            io_writers: 1,
            buffer_blocks: 10,
            block_size_override: Some(5),
            direct_io: false,
            max_buffer_memory: 100_000_000,
        };

        rt.block_on(async {
            let handle = store_and_drop_matrix_streaming_optimized(
                matrix,
                temp_dir.path(),
                "test_matrix",
                config,
            );

            handle.await.unwrap();
        });

        // Verify files were created
        let entries: Vec<_> = fs::read_dir(temp_dir.path()).unwrap().collect();
        assert!(!entries.is_empty(), "No matrix files were created");

        // Check that matrix files exist
        let matrix_files: Vec<_> = entries
            .into_iter()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_name().to_string_lossy().contains("test_matrix"))
            .collect();

        assert!(!matrix_files.is_empty(), "No test_matrix files found");
    }

    #[test]
    fn test_performance_comparison() {
        let rt = Runtime::new().unwrap();
        let temp_dir1 = TempDir::new().unwrap();
        let temp_dir2 = TempDir::new().unwrap();

        // Create identical matrices for comparison
        let matrix1 = create_test_matrix(50, 50);
        let matrix2 = create_test_matrix(50, 50);

        let optimized_config = StorageConfig {
            serialization_workers: 4,
            io_writers: 2,
            buffer_blocks: 50,
            block_size_override: Some(25),
            direct_io: false,
            max_buffer_memory: 500_000_000,
        };

        rt.block_on(async {
            // Test original implementation
            let start1 = Instant::now();
            let handle1 = crate::storage::store_and_drop_matrix_streaming(
                matrix1,
                temp_dir1.path(),
                "original_matrix",
            );
            handle1.await.unwrap();
            let original_time = start1.elapsed();

            // Test optimized implementation
            let start2 = Instant::now();
            let handle2 = store_and_drop_matrix_streaming_optimized(
                matrix2,
                temp_dir2.path(),
                "optimized_matrix",
                optimized_config,
            );
            handle2.await.unwrap();
            let optimized_time = start2.elapsed();

            println!("Original implementation: {:?}", original_time);
            println!("Optimized implementation: {:?}", optimized_time);

            // Verify both created files
            let original_files = fs::read_dir(temp_dir1.path()).unwrap().count();
            let optimized_files = fs::read_dir(temp_dir2.path()).unwrap().count();

            assert!(original_files > 0, "Original implementation created no files");
            assert!(optimized_files > 0, "Optimized implementation created no files");
            assert_eq!(original_files, optimized_files, "Different number of files created");
        });
    }

    #[test]
    fn test_large_buffer_memory_management() {
        let rt = Runtime::new().unwrap();
        let temp_dir = TempDir::new().unwrap();
        let matrix = create_test_matrix(20, 20);

        let config = StorageConfig {
            serialization_workers: 2,
            io_writers: 1,
            buffer_blocks: 5, // Small buffer to test memory management
            block_size_override: Some(1000),
            direct_io: false,
            max_buffer_memory: 10_000_000, // Small memory limit
        };

        rt.block_on(async {
            let handle = store_and_drop_matrix_streaming_optimized(
                matrix,
                temp_dir.path(),
                "memory_test",
                config,
            );

            // Should complete without panicking despite memory constraints
            handle.await.unwrap();
        });

        // Verify files were still created despite memory pressure
        let files = fs::read_dir(temp_dir.path()).unwrap().count();
        assert!(files > 0, "No files created under memory pressure");
    }

    #[test]
    fn test_multiple_io_writers() {
        init_tracing();
        let rt = Runtime::new().unwrap();
        let temp_dir = TempDir::new().unwrap();
        let matrix = create_test_matrix(100, 100);

        let config = StorageConfig::default();
        rt.block_on(async {
            let handle = store_and_drop_matrix_streaming_optimized(
                matrix,
                temp_dir.path(),
                "multi_io_test",
                config,
            );

            handle.await.unwrap();
        });

        let files = fs::read_dir(temp_dir.path()).unwrap().count();
        assert!(files > 0, "Multiple I/O writers failed to create files");
    }
}
