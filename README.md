# diamond-io

Implementation of [Diamond iO](https://eprint.iacr.org/2025/236)

### Note

We currently support two different matrix implementations:
1. **In-memory** (default): Uses memory for all matrix storage.
2. **Disk-backed** (enable with `--features disk`): Uses the `mmap()` syscall to store matrices on disk.

#### Test iO (without `test` feature)

This disables helper logic and fields used only for testing, which are not required for iO security.

- **Dummy parameters**  
```bash
cargo test -r --test test_io_dummy_param --no-default-features -- --ignored --nocapture
```

- **Real parameters** (tests are ignored by default)  
```bash
cargo test -r --test test_io_real_param --no-default-features -- --ignored --nocapture
```

- **With memory profiler**  
```bash
uv run memory_profile.py cargo test -r --test test_io_dummy_param --no-default-features
```