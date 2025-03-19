## diamond-io

Implementation of [Diamond iO](https://eprint.iacr.org/2025/236)

#### Full Test (with test feature)

```
cargo test -r  --features="parallel"
```

#### Test iO (without test feature) 
this will remove helper logic + helper fields for test which is not require for iO security.

```
cargo test --test test_io --no-default-features --features parallel
```