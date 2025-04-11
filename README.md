## diamond-io

Implementation of [Diamond iO](https://eprint.iacr.org/2025/236)

#### Full Test (with test feature)

```
cargo test -r
```

#### Test iO (without test feature) 

this will remove helper logic + helper fields for test which is not require for iO security.


- dummy params
```
cargo test -r --test test_io_dummy_param test_io_just_mul_enc_and_bit --no-default-features -- --nocapture
```

- real params (by default ignored)
```
cargo test -r --test test_io_real_param --no-default-features -- --ignored --nocapture
```

- with memory profiler 
```
uv run memory_profile.py cargo test -r --test test_io_dummy_param --no-default-features
```