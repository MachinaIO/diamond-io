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
            "323778148704285877904461387615990672855113714166872436505644678160965660026126816597",
            1,
            1,
            1,
            12.057,
            542800000000000000000000000.0,
            "tests/io_middle_param",
        )
        .await;
    }
}
