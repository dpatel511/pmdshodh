//! Summary JSON: arbitrary bytes through `report::parse` and `render`. Expect an error naming the
//! file, or HTML whose numeric attributes (SVG coordinates, plotted values, thresholds) are all
//! finite numbers; never a panic. Text may repeat user input (file names, skip reasons), so only
//! numeric attributes are checked.
#![no_main]

use libfuzzer_sys::fuzz_target;
use pmdshodh::report;

const NUMERIC: [&str; 12] = [
    "x", "y", "x1", "x2", "y1", "y2", "cx", "cy", "width", "height", "points", "data-threshold",
];

/// Values of the numeric attributes in `html`.
fn numeric_values(html: &str) -> Vec<&str> {
    let mut out = Vec::new();
    for name in NUMERIC {
        let key = format!(" {name}=\"");
        let mut rest = html;
        while let Some(i) = rest.find(&key) {
            let after = &rest[i + key.len()..];
            let end = after.find('"').unwrap();
            out.extend(after[..end].split([' ', ',']).filter(|v| !v.is_empty()));
            rest = &after[end..];
        }
    }
    // data-values: numbers, "n/a", or "length:count" pairs.
    let mut rest = html;
    while let Some(i) = rest.find(" data-values=\"") {
        let after = &rest[i + 14..];
        let end = after.find('"').unwrap();
        out.extend(
            after[..end]
                .split([' ', ':'])
                .filter(|v| !v.is_empty() && *v != "n/a"),
        );
        rest = &after[end..];
    }
    out
}

fuzz_target!(|data: &[u8]| {
    let Ok(json) = serde_json::from_slice::<serde_json::Value>(data) else {
        return;
    };
    match report::parse("f.json", &json) {
        Err(e) => assert!(e.starts_with("f.json: "), "{e}"),
        Ok(mut sample) => {
            sample.label = "s".into();
            let html = report::render(&[sample]);
            for v in numeric_values(&html) {
                let n: f64 = v.parse().unwrap_or_else(|_| panic!("not a number: {v:?}"));
                assert!(n.is_finite(), "{v}");
            }
        }
    }
});
