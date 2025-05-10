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
            17,
            "182270893731917660364375185298660602983033397643884786244994464562915433787344833915025064536484522",
            1,
            1,
            1,
            12.919,
            108910000000000000000000000000000000.0,
            "tests/io_real_param",
        )
        .await;
    }
}
