#[cfg(test)]
mod test {
    use diamond_io::test_utils::test_io_common;

    #[tokio::test]
    #[ignore]
    async fn test_io_just_mul_enc_and_and_bit() {
        test_io_common(4, 2, 17, 10, "1", 3, 4, 1, 0.0, 0.0, "tests/io_dummy_param").await;
    }
}
