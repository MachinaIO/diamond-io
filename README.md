# diamond-io

Implementation of [Diamond iO](https://eprint.iacr.org/2025/236), a straightforward construction of indistinguishability obfuscation (iO).

## Installation

### Prerequisites
- [Rust](https://www.rust-lang.org/tools/install) 1.87 nightly
- [OpenFHE](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/installation.html) (System install required in `/usr/local/lib`), make sure to install our [fork](https://github.com/MachinaIO/openfhe-development/tree/feat/improve_determinant) in `/feat/improve_determinant` branch

## Experiments

You can run evaluation experiments with [e2e](./e2e/) parameters with [dio](./dio/) cli tooling.

### Building

After installing the prerequisites, you can build the project using:

```bash
cargo build --release
```

## Test iO (End-To-End)

- **Dummy parameters**  
Fastest way to check if the end-to-end process works with insecure parameters:
```bash
cargo test -r --test test_io_dummy_param --no-default-features -- --nocapture
```

- **Real parameters** 
Warning: You need sufficient RAM.
```bash
cargo test -r --test test_io_real_param --no-default-features -- --ignored --nocapture
```

- **With memory profiler**  
```bash
uv run memory_profile.py cargo test -r --test test_io_dummy_param --no-default-features
```

## Note

We currently support two different matrix implementations:
1. **In-memory** (default): Uses memory for all matrix storage.
2. **Disk-backed** (enable with `--features disk`): Uses the `mmap()` syscall to store matrices on disk.


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
Among the parameters, `crt_depth` denotes the minimum number of moduli satisfying correctness and security, and `d`, `hardcoded_key_sigma`, `p_sigma`, and `switched_modulus` can be used for `ObfuscationParams`.


## Acknowledgments

*We would like to sincerely thank the developers of [OpenFHE](https://github.com/openfheorg/openfhe-development) and [openfhe-rs](https://github.com/fairmath/openfhe-rs), open-source lattice and FHE libraries, whose optimized implementations of trapdoor sampling, RLWE primitives, and Rust bindings played a crucial role in helping us implement Diamond iO. We are also grateful to Prof. Yuriy Polyakov for his valuable advice on preimage sampling and his insightful feedback on optimizing our implementation. We greatefully acknowledge [Community Privacy Residency (2025)](https://community-privacy.github.io/partners/), in which our earliest implementation was developed. Any remaining errors are entirely our own responsibility.*