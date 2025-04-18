#[cfg(test)]
mod test {
    use diamond_io::test_utils::test_io_common;

    #[tokio::test]
    async fn test_io_just_mul_enc_and_bit_middle_params() {
        test_io_common(
            4096,
            5,
            51,
            20,
            "431359146674410224742050828377557384853608161148858662676258371928064",
            1,
            1,
            1,
            15.82764,
            15.82764,
            15.82764,
            "tests/io_middle_param",
        )
        .await;
    }
}
