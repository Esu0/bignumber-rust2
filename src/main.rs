mod bignum;

use rand::Rng;
use std::time::Instant;

fn main() {
    time();
}

fn time() {
    let n = 1 << 22;
    let mut rng = rand::thread_rng();
    let a: Vec<i64> = (0..n).map(|_| rng.gen_range(0..1000000000)).collect();
    let b: Vec<i64> = (0..n).map(|_| rng.gen_range(0..1000000000)).collect();
    let s = Instant::now();
    let v = bignum::number_theoretic::convolution::convolution_level2(a.into_iter(), b.into_iter(), n, n);
    println!("elapsed: {} ms", s.elapsed().as_micros() as f32 / 1000.);
}