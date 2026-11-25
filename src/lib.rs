//! pmdshodh: per-read post-mortem damage (PMD) scoring and filtering of ancient DNA alignments.
//!
//! pmdshodh is a command-line tool; install it with `cargo install pmdshodh` and see the
//! [README](https://github.com/dpatel511/pmdshodh#readme). Its stable interfaces are the command
//! line, its outputs and the JSON run summary (`format_version`). This library holds the
//! command's internals; its API may change in any release, including patch releases.

pub mod bgzf;
pub mod cram;
pub mod damage;
pub mod model;
pub mod reads;
pub mod reference;
pub mod report;
pub mod sites;
pub mod summary;
