// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use std::f64::consts::PI;
const RES: f64 = 1e-10;
const MINNUMLOOP: i32 = 3;

pub fn poisson_pval(x:i32, lambda:f64) -> f64 {
    let mut sum = 0.0;
    let mut i = x;
    loop {
        let term = lambda.powf(i as f64) / fast_factorial(i as f64);
        sum += term;
        if term < RES && i - x >=MINNUMLOOP { break}
        i += 1;
    }
    (-lambda).exp()*sum
}

pub fn factorial(x:u16) -> f64 {
    match x {
        0 => 1.0,
        1 => 1.0,
        _ => x as f64*factorial(x-1),
    }
}

pub fn fast_factorial(x:f64) -> f64 {
    let a: f64 = (2.0*x+1.0/3.0)*PI;
    a.sqrt()*x.powf(x)*(-x).exp()
}
