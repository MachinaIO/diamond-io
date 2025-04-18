# Set environment variable
export CARGO_MAKE_EXTEND_WORKSPACE_MAKEFILE := "true"

# Run rustfmt to check the code formatting without making changes
format:
    cargo +nightly fmt -- --check

# Clean up the project by removing the target directory
clean:
    cargo clean

# Run clippy to catch common mistakes and improve your Rust code
clippy:
    RUSTFLAGS="-A unused" cargo +nightly clippy --all-targets --all-features -- -Dwarnings

# Generate documentation for the project
docs:
    cargo doc --no-deps

# Execute all unit tests in the workspace
test:
   cargo test -r

test-io:
   cargo test -r --test test_io_dummy_param --no-default-features

e2e:
    dio run -c e2e/dio-config.dummy.toml --add-num 0 --mul-num 1
    dio run -c e2e/dio-config.dummy.toml --add-num 1 --mul-num 1
    dio run -c e2e/dio-config.dummy.toml --add-num 1 --mul-num 2
    dio run -c e2e/dio-config.dummy.toml --add-num 1 --mul-num 3
    dio run -c e2e/dio-config.dummy.toml --add-num 2 --mul-num 1
    dio run -c e2e/dio-config.dummy.toml --add-num 3 --mul-num 1
    dio run -c e2e/dio-config.dummy.toml --add-num 4 --mul-num 1
    dio run -c e2e/dio-config.dummy.toml --add-num 5 --mul-num 1
    dio run -c e2e/dio-config.dummy-width.toml --add-num 0 --mul-num 1
    dio run -c e2e/dio-config.dummy-width.toml --add-num 1 --mul-num 1
    dio run -c e2e/dio-config.dummy-width.toml --add-num 1 --mul-num 2
    dio run -c e2e/dio-config.dummy-width.toml --add-num 1 --mul-num 3
    dio run -c e2e/dio-config.dummy-width.toml --add-num 2 --mul-num 1
    dio run -c e2e/dio-config.dummy-width.toml --add-num 3 --mul-num 1
    dio run -c e2e/dio-config.dummy-width.toml --add-num 4 --mul-num 1
    dio run -c e2e/dio-config.dummy-width.toml --add-num 5 --mul-num 1

# Run the entire CI pipeline including format, clippy, docs, and test checks
ci: format clippy docs test test-io
    @echo "CI flow completed"