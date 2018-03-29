// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use std::f64::consts::PI;
use special::Gamma;

pub fn ln_stirling_gosper(x: i32) -> f64 {
    if x==0 {return 0.0f64;}
    // Gosper's Stierling approximation for better accuracy
    let n = x as f64;
    n * n.ln() - n + 0.5 * ((2.0*n+1.0/3.0)* PI).ln()
}

pub fn ln_factorial(x: f64) -> f64 {
    (x + 1.0).ln_gamma().0
}


pub fn loglikelihood(exp: f64, obs: i32) -> f64 {
    let exp1 = if exp < 0.1f64 {0.1f64} else {exp};
    -exp1 + obs as f64* exp1.ln() - ln_factorial(obs as f64)
}


pub fn binomial_coefficient_log(n: f64, k: f64) -> f64 {
        ln_factorial(n) - ln_factorial(k) - ln_factorial(n-k)
}

pub fn negative_binomial_loglikelihood(mu: f64, alpha: f64, obs: i32) -> f64 {
    assert!(alpha > 1.0);

    let p = 1.0 - 1.0/alpha;
    let r = mu / (alpha - 1.0); 
    let k = obs as f64;

    binomial_coefficient_log(k+ r - 1.0, k) + r * (1.0/alpha).ln() + k * p.ln()

}


#[cfg(test)]
mod tests {
    use assert;
    use super::*;

    #[test]
    fn poisson() {
        assert::close(loglikelihood(1.0, 10), -16.10441, 1e-5);
    }

    #[test]
    fn neg_bin() {
        assert::close(negative_binomial_loglikelihood(1.0, 1.01, 10), -15.76, 1e-1);
    }
}
