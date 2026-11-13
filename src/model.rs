//! Per-site PMD score terms, matching PMDtools 0.60.
//!
//! Built from the paper (Skoglund et al. 2014, PNAS, doi:10.1073/pnas.1318934111) and from
//! black-box tests of PMDtools 0.60 (docs/blackbox.md section 1), not from PMDtools' source.
//! A read's PMD score is the sum of the terms of its informative sites.
//! The readable reference version is `tests/oracle/pmd_reference.py`.

/// Damage assumed by the null model. The paper uses 0; PMDtools 0.60 uses 0.001 and has no
/// option to change it (docs/blackbox.md section 2).
pub const NULL_BACKGROUND: f64 = 0.001;

/// Model parameters, named after the PMDtools options that set them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Params {
    /// `--PMDpparam`: geometric decay of damage away from the read end.
    pub p: f64,
    /// `--PMDconstant`: damage floor added at every position.
    pub c: f64,
    /// `--polymorphism_ancient`: polymorphism rate in the ancient (PMD) model.
    pub pi_ancient: f64,
    /// `--polymorphism_contamination`: polymorphism rate in the null model.
    pub pi_contamination: f64,
}

impl Default for Params {
    /// PMDtools 0.60 defaults.
    fn default() -> Self {
        Params {
            p: 0.3,
            c: 0.01,
            pi_ancient: 0.001,
            pi_contamination: 0.001,
        }
    }
}

/// Probability that a sequencing error gives the one modelled mismatch: 10^(-Q/10) / 3.
/// No cap on Q, as in PMDtools.
pub fn error_rate(q: u8) -> f64 {
    10f64.powf(-f64::from(q) / 10.0) / 3.0
}

/// Damage probability at distance `z` (1 = terminal base) from the relevant end (paper Eq. 3).
pub fn damage(z: u32, p: f64, c: f64) -> f64 {
    (1.0 - p).powf(f64::from(z) - 1.0) * p + c
}

/// Probability of a match at a site with damage `d` (paper Eq. 1).
fn p_match(d: f64, pi: f64, e: f64) -> f64 {
    (1.0 - pi) * (1.0 - e) * (1.0 - d) + (1.0 - pi) * e * d + pi * e * (1.0 - d)
}

/// Match probability with damage `a` and `b` from the two ends (paper SI Eq. S1).
fn p_match_two_ends(a: f64, b: f64, pi: f64, e: f64) -> f64 {
    (1.0 - pi) * (1.0 - e) * (1.0 - a) * (1.0 - b)
        + (1.0 - pi) * e * a * (1.0 - b)
        + (1.0 - pi) * e * b * (1.0 - a)
        + pi * e * (1.0 - a) * (1.0 - b)
}

/// Log-likelihood ratio of one observation, given each model's match probability.
fn log_ratio(is_match: bool, p_ancient: f64, p_null: f64) -> f64 {
    if is_match {
        (p_ancient / p_null).ln()
    } else {
        ((1.0 - p_ancient) / (1.0 - p_null)).ln()
    }
}

/// Double-stranded score term of one informative site at distance `z` with base quality `q`.
/// `is_match` is false for a C→T (or G→A) mismatch.
pub fn site_term_ds(is_match: bool, z: u32, q: u8, params: &Params) -> f64 {
    let e = error_rate(q);
    log_ratio(
        is_match,
        p_match(damage(z, params.p, params.c), params.pi_ancient, e),
        p_match(NULL_BACKGROUND, params.pi_contamination, e),
    )
}

/// Single-stranded (`--ss`) score term of one C site, `z` from the 5′ end and `y` from the 3′ end,
/// as PMDtools 0.60 computes it: the double-stranded term plus the S1 term, so each site counts
/// twice (docs/blackbox.md section 1.4). Used only with `--pmdtools-compat`.
pub fn site_term_ss(is_match: bool, z: u32, y: u32, q: u8, params: &Params) -> f64 {
    site_term_ds(is_match, z, q, params) + site_term_s1(is_match, z, y, q, params)
}

/// Single-stranded score term of one C site from paper SI Eq. S1 alone (the default `--ss`
/// score), with the same null-model background damage as the double-stranded term.
pub fn site_term_s1(is_match: bool, z: u32, y: u32, q: u8, params: &Params) -> f64 {
    let e = error_rate(q);
    log_ratio(
        is_match,
        p_match_two_ends(
            damage(z, params.p, params.c),
            damage(y, params.p, params.c),
            params.pi_ancient,
            e,
        ),
        p_match_two_ends(NULL_BACKGROUND, NULL_BACKGROUND, params.pi_contamination, e),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn error_rate_is_a_third_of_phred_probability() {
        assert!(close(error_rate(0), 1.0 / 3.0));
        assert!(close(error_rate(10), 0.1 / 3.0));
        assert!(close(error_rate(40), 1e-4 / 3.0));
    }

    #[test]
    fn damage_follows_paper_eq3_with_defaults() {
        assert!(close(damage(1, 0.3, 0.01), 0.31));
        assert!(close(damage(2, 0.3, 0.01), 0.22));
        assert!(close(damage(3, 0.3, 0.01), 0.157));
    }

    // Values printed by PMDtools 0.60 --printDS for one-site reads (docs/blackbox.md, probe1/probe6).
    #[test]
    fn ds_terms_match_pmdtools_at_q40() {
        let p = Params::default();
        assert!(close(site_term_ds(false, 1, 40, &p), 5.029717285061266));
        assert!(close(site_term_ds(true, 1, 40, &p), -0.37004823819308086));
        assert!(close(site_term_ds(false, 30, 40, &p), 1.6917038398804485));
    }

    #[test]
    fn ss_terms_match_pmdtools_at_q40() {
        let p = Params::default();
        assert!(close(site_term_ss(false, 1, 30, 40, &p), 9.681909363978177));
        assert!(close(site_term_ss(false, 30, 1, 40, &p), 6.343895918797358));
        assert!(close(
            site_term_ss(true, 1, 30, 40, &p),
            -0.7491557655766737
        ));
    }

    // Departure from PMDtools (docs/design.md §5): the default --ss term is Eq. S1 alone, one
    // term per site and symmetric in the two ends. Values computed independently from the
    // Eq. S1 formula (docs/blackbox.md 1.4) with the null background 0.001.
    #[test]
    fn s1_term_counts_each_site_once() {
        let p = Params::default();
        assert!(close(site_term_s1(false, 1, 30, 40, &p), 4.65219207891691));
        assert!(close(site_term_s1(false, 30, 1, 40, &p), 4.65219207891691));
        assert!(close(
            site_term_s1(true, 1, 30, 40, &p),
            -0.37910752738359293
        ));
        // PMDtools' value is the ds term plus this one.
        assert!(close(
            site_term_ss(false, 1, 30, 40, &p),
            site_term_ds(false, 1, 40, &p) + site_term_s1(false, 1, 30, 40, &p)
        ));
    }
}
