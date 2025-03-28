use rodio::{Decoder, OutputStream, OutputStreamBuilder, Sink, Source};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub struct Player {
    stream_handle: OutputStream,
    sink: Sink,
}

impl Player {
    pub fn new() -> Self {
        let stream_handle = OutputStreamBuilder::open_default_stream().unwrap();
        let sink = rodio::Sink::connect_new(&stream_handle.mixer());
        Player { stream_handle, sink }
    }

    pub fn play_music<P: AsRef<Path>>(&self, file: P) {
        self.stop_music();
        let file = File::open(file).unwrap();
        let source = rodio::Decoder::try_from(file).unwrap().repeat_infinite();
        self.sink.append(source);
        self.sink.play();
    }

    pub fn stop_music(&self) {
        self.sink.stop();
        self.sink.clear();
    }
}
