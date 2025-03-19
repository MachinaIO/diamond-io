## diamond-io

Implementation of [Diamond iO](https://eprint.iacr.org/2025/236)

#### Test (with test feature)

```
cargo test  --features="parallel"
```

#### Test (without test feature) 
- io test

```
cargo test --test test_io --no-default-features --features parallel
```



run io in parallel
```
cargo test -r --package diamond-io --lib --features="parallel" -- io::test::test_io_just_mul_enc_and_bit --exact --show-output
```