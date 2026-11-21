//! Summary JSON: arbitrary bytes through `report::parse` and `render`. Expect an error naming the
//! file, or HTML with no NaN or infinity in it; never a panic.
#![no_main]

use libfuzzer_sys::fuzz_target;
use pmdshodh::report;

fuzz_target!(|data: &[u8]| {
    let Ok(json) = serde_json::from_slice::<serde_json::Value>(data) else {
        return;
    };
    match report::parse("f.json", &json) {
        Err(e) => assert!(e.starts_with("f.json: "), "{e}"),
        Ok(mut sample) => {
            sample.label = "s".into();
            let html = report::render(&[sample]);
            assert!(!html.contains("NaN"));
            for bad in ["inf\"", "inf,", "inf ", "inf<"] {
                assert!(!html.contains(bad), "{bad}");
            }
        }
    }
});
