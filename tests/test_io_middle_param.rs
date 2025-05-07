#[cfg(test)]
mod test {
    use diamond_io::test_utils::test_io_common;

    #[tokio::test]
    async fn test_io_just_mul_enc_and_bit_middle_params() {
        test_io_common(
            4096,
            5,
            51,
            17,
            "9202328462379048914265825575350919550703868407672378932957413585715200",
            1,
            1,
            1,
            15.828,
            372250000000.0,
            "tests/io_middle_param",
        )
        .await;
    }
}
