# diamond-io

Implementation of [Diamond iO](https://eprint.iacr.org/2025/236)

## Note

We currently support two different matrix implementations:
1. **In-memory** (default): Uses memory for all matrix storage.
2. **Disk-backed** (enable with `--features disk`): Uses the `mmap()` syscall to store matrices on disk.

## Test iO (without `test` feature)

This disables helper logic and fields used only for testing, which are not required for iO security.

- **Dummy parameters**  
```bash
cargo test -r --test test_io_dummy_param --no-default-features -- --nocapture
```

- **Real parameters** (tests are ignored by default)  
```bash
cargo test -r --test test_io_real_param --no-default-features -- --ignored --nocapture
```

- **With memory profiler**  
```bash
uv run memory_profile.py cargo test -r --test test_io_dummy_param --no-default-features
```

## Simulate Parameters
Our simulator only targets circuits used for our benchmarks.

1. Make sure to install [`dio`](/dio/) binary before
2. Change the following values hardcoded in `simulator/main.py` after the line `if __name__ == "__main__":`:
    - `secpar`: the minimum security parameter you want to guarantee.
    - `log2_n`: a log2 value of the ring dimension.
    - `max_d`: the maximum value of the number of the secret polynomials denoted by `d`.
    - `min_base_bits`: the minimum value of the base bits for decomposition denoted by `base_bits`.
    - `max_base_bits`: the maximum value of `base_bits`.
    - `crt_bits`: the bits of each moduli of CRT.
    - `max_crt_depth`: the maximum number of moduli.
    - `input_size`: the evaluator's input bit size.
    - `input_width`: the number of bits inserted at each diamond. The larger value of `input_width` increase the number of preimages but decrease the required modulus size.
    - `add_num`: the number of addition gates for the evaluator's input bits.
    - `mul_num`: the number of multiplication gates for the evaluator's input bits.
3. Install sagemath if you have not installed it. Ref: https://doc.sagemath.org/html/en/installation/conda.html
4. Run `sage main.py` under the `simulator` directory.

If the script is completed without any error, the found parameters are added to the last line in `simulator/params.log`. 
Among the parameters, `crt_depth` denotes the minimum number of moduli satisfying correctness and security, and `d`, `encoding_sigma`, `hardcoded_key_sigma`, `p_sigma`, and `switched_modulus` can be used for `ObfuscationParams`.
