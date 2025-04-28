# dio

diamond-io CLI implementation

### install
```
cargo install --path .
```
disk feature enable
```
cargo install --features disk --path .
```
debug feature enable
```
cargo install --features debug --path .
```

### Run

Execute the end-to-end workflow from obfuscation to evaluation:

```
dio run-bench \
  -c {CONFIG-TOML-PATH} \
  -o {OBFUSCATION-DIRECTORY-PATH} \
  --add-num {ADD-GATE-NUMBER} \
  --mul-num {MUL-GATE-NUMBER}
```

### Build Circuit

Build the final circuit from the specified gate numbers and configuration, and count its gate types:

```
dio build-circuit \
  -c {CONFIG-TOML-PATH} \
  --add-num {ADD-GATE-NUMBER} \
  --mul-num {MUL-GATE-NUMBER}
```