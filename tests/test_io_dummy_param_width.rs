#[cfg(test)]
mod test {
    use diamond_io::test_utils::test_io_common;

    #[tokio::test]
    async fn test_io_just_mul_enc_and_and_bit_width() {
        test_io_common(4, 2, 17, 10, "1", 3, 4, 2, 0.0, 0.0, 0.0, "tests/io_dummy_param_width")
            .await;
    }
}
