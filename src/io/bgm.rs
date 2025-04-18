use rodio::{OutputStreamBuilder, Sink, Source};
use std::{fs::File, path::Path};

pub struct Player {
    sink: Sink,
}

impl Default for Player {
    fn default() -> Self {
        Self::new()
    }
}

impl Player {
    pub fn new() -> Self {
        let stream_handle = OutputStreamBuilder::open_default_stream().unwrap();
        let sink = rodio::Sink::connect_new(stream_handle.mixer());
        Player { sink }
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
