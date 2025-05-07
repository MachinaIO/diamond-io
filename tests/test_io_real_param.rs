#[cfg(test)]
mod test {
    use diamond_io::test_utils::test_io_common;

    #[tokio::test]
    #[ignore]
    async fn test_io_just_mul_enc_and_bit_real_params() {
        test_io_common(
            8192,
            6,
            51,
            17,
            "7286787346663540544401205006872800924021905474424898872905693286510586401192117719085",
            1,
            1,
            1,
            12.057,
            40616000000000000000.0,
            "tests/io_real_param",
        )
        .await;
    }
}
