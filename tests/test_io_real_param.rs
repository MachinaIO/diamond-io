#[cfg(test)]
mod test {
    use diamond_io::test_utils::test_io_common;

    #[tokio::test]
    #[ignore]
    async fn test_io_just_mul_enc_and_bit_real_params() {
        test_io_common(
            8192,
            7,
            51,
            20,
            "546812681195752981093125556779405341132604288058152353837071503782819436677001997601680811337187328",
            1,
            1,
            1,
            12.91885,
            108910396484176728921799104269415406.33545,
            12.91885,
            "tests/io_real_param",
        )
        .await;
    }
}
