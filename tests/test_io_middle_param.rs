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
            "431359140246658059678711139224465721413543900671837184566486422781952",
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
