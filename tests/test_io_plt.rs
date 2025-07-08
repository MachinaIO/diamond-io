#[cfg(test)]
mod test {
    use diamond_io::test_utils::test_io_plt;

    #[tokio::test]
    async fn test_io_plt_dummy() {
        test_io_plt(4, 2, 17, 10, "1", 3, 3, 1, 0.0, 0.0, "tests/io_plt").await;
    }
}
