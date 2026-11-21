//! BGZF: arbitrary bytes through the reader (an error or data, never a panic or hang), and any
//! data written and read back unchanged, with and without the thread pool.
#![no_main]

use std::io::{Read, Write};

use libfuzzer_sys::fuzz_target;
use pmdshodh::bgzf::{Reader, Writer};

fuzz_target!(|data: &[u8]| {
    let parallel = data.first().is_some_and(|b| b & 1 == 1);
    let mut out = Vec::new();
    let _ = Reader::new(data, parallel).read_to_end(&mut out);

    let mut compressed = Vec::new();
    {
        let mut w = Writer::new(&mut compressed, parallel);
        w.write_all(data).unwrap();
        w.finish().unwrap();
    }
    let mut back = Vec::new();
    Reader::new(&compressed[..], parallel)
        .read_to_end(&mut back)
        .unwrap();
    assert_eq!(back, data);
});
