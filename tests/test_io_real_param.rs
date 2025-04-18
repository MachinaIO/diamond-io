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
            20,
            "242833611528216133864932738352844082358996736827870043467279656893386864455514587136",
            1,
            1,
            1,
            12.05698,
            40615715852990820734.97011,
            12.05698,
            "tests/io_real_param",
        )
        .await;
    }
}
