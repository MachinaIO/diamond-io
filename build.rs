fn main() {
    println!("cargo:rerun-if-changed=src/main.rs");

    // linking openFHE
    println!("cargo:rustc-link-arg=-L/usr/local/lib");
    println!("cargo:rustc-link-arg=-lOPENFHEpke");
    println!("cargo:rustc-link-arg=-lOPENFHEbinfhe");
    println!("cargo:rustc-link-arg=-lOPENFHEcore");

    // linking OpenMP
    println!("cargo:rustc-link-arg=-fopenmp");

    // necessary to avoid LD_LIBRARY_PATH
    println!("cargo:rustc-link-arg=-Wl,-rpath,/usr/local/lib");

    // Only execute the following code when the "bgm" feature is enabled
    #[cfg(feature = "bgm")]
    {
        use std::fs;
        use std::path::Path;
        // Create bgm directory if it doesn't exist
        let bgm_dir = Path::new("bgm");
        if !bgm_dir.exists() {
            println!("Creating bgm directory");
            fs::create_dir(bgm_dir).expect("Failed to create bgm directory");
        }

        // Download bgm/obf_bgm1.mp3 if it doesn't exist
        let bgm_files = vec![
            ("obf_bgm1.mp3", "https://bgmer.net/wp-content/uploads/2022/03/233_BPM163.mp3"),
            ("obf_bgm2.mp3", "https://bgmer.net/wp-content/uploads/2024/02/420_BPM108.mp3"),
            ("obf_bgm3.mp3", "https://bgmer.net/wp-content/uploads/2021/12/251_BPM150.mp3"),
            ("obf_bgm4.mp3", "https://bgmer.net/wp-content/uploads/2021/09/157_BPM175.mp3"),
            ("obf_bgm5.mp3", "https://bgmer.net/wp-content/uploads/2021/05/075_BPM140.mp3"),
        ];
        for (file_name, url) in bgm_files {
            let bgm_file = bgm_dir.join(file_name);
            if !bgm_file.exists() {
                println!("Downloading {}", file_name);
                let response = reqwest::blocking::get(url).expect("Failed to download file");
                let content = response.bytes().expect("Failed to read response");
                fs::write(bgm_file, content).expect("Failed to write file");
                println!("Downloaded {}", file_name);
            }
        }
    }
}
