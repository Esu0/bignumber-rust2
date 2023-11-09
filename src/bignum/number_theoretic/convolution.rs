#![allow(dead_code)]
use super::{NttTable, ParallelConvoluter};
use itertools::Itertools;
use std::thread;

struct Tables(
    NttTable<2013265921, 31>,
    NttTable<1811939329, 13>,
    NttTable<469762049, 3>,
);

impl Tables {
    const MOD: (i64, i64, i64) = (2013265921, 1811939329, 469762049);
    const fn new() -> Self {
        Self(NttTable::new(), NttTable::new(), NttTable::new())
    }
}

impl<'a> ParallelConvoluter for &'a Tables {
    type Ret = (Vec<i64>, Vec<i64>, Vec<i64>);
    fn convolution_parallel(
        self,
        a: impl Iterator<Item = i64> + Send,
        b: impl Iterator<Item = i64>,
        size_a: usize,
        size_b: usize,
    ) -> Self::Ret {
        let s = size_a + size_b - 1;
        let mut va: (Vec<i64>, Vec<i64>, Vec<i64>) = a
            .take(s)
            .map(|n| {
                (
                    n.rem_euclid(Tables::MOD.0),
                    n.rem_euclid(Tables::MOD.1),
                    n.rem_euclid(Tables::MOD.2),
                )
            })
            .multiunzip();
        let mut vb: (Vec<i64>, Vec<i64>, Vec<i64>) = b
            .take(s)
            .map(|n| {
                (
                    n.rem_euclid(Tables::MOD.0),
                    n.rem_euclid(Tables::MOD.1),
                    n.rem_euclid(Tables::MOD.2),
                )
            })
            .multiunzip();
        let ex = s.next_power_of_two();
        thread::scope(|s| {
            s.spawn(|| {
                thread::scope(|s2| {
                    s2.spawn(|| {
                        va.1.resize(ex, 0);
                        self.1.ntt_mut_slc_unchecked(&mut va.1);
                    });
                    vb.1.resize(ex, 0);
                    self.1.ntt_mut_slc_unchecked(&mut vb.1);
                });
                // self.1.ntt_mut_slc_unchecked(&mut va.1);
                va.1.iter_mut()
                    .zip(vb.1)
                    .for_each(|(an, bn)| *an = *an * bn % Tables::MOD.1);
                let invsize = super::minv(ex as u32, Tables::MOD.1 as u32) as i64;
                va.1.iter_mut()
                    .for_each(|an| *an = *an * invsize % Tables::MOD.1);
                self.1.intt_mut_slc_unchecked(&mut va.1);
            });
            s.spawn(|| {
                thread::scope(|s2| {
                    s2.spawn(|| {
                        va.2.resize(ex, 0);
                        self.2.ntt_mut_slc_unchecked(&mut va.2);
                    });
                    vb.2.resize(ex, 0);
                    self.2.ntt_mut_slc_unchecked(&mut vb.2);
                });
                va.2.iter_mut()
                    .zip(vb.2)
                    .for_each(|(an, bn)| *an = *an * bn % Tables::MOD.2);
                let invsize = super::minv(ex as u32, Tables::MOD.2 as u32) as i64;
                va.2.iter_mut()
                    .for_each(|an| *an = *an * invsize % Tables::MOD.2);
                self.2.intt_mut_slc_unchecked(&mut va.2);
            });
            thread::scope(|s2| {
                s2.spawn(|| {
                    va.0.resize(ex, 0);
                    self.0.ntt_mut_slc_unchecked(&mut va.0);
                });
                vb.0.resize(ex, 0);
                self.0.ntt_mut_slc_unchecked(&mut vb.0);
            });
            va.0.iter_mut()
                .zip(vb.0)
                .for_each(|(an, bn)| *an = *an * bn % Tables::MOD.0);
            let invsize = super::minv(ex as u32, Tables::MOD.0 as u32) as i64;
            va.0.iter_mut()
                .for_each(|an| *an = *an * invsize % Tables::MOD.0);
            self.0.intt_mut_slc_unchecked(&mut va.0);
        });
        va
    }
}
//const NTTPRIMROOT: (i64, i64, i64) = (31, 13, 3);
const NTT_TABLE: Tables = Tables::new();

pub fn convolution_mod1_mut_unchecked(a: &mut [i64], b: &mut [i64]) {
    NTT_TABLE.0.ntt_mut_slc_unchecked(a);
    NTT_TABLE.0.ntt_mut_slc_unchecked(b);
    a.iter_mut()
        .zip(b.iter())
        .for_each(|(an, bn)| *an = *an * bn % Tables::MOD.0);
    NTT_TABLE.0.intt_mut_slc_unchecked(a);
    let invsize = super::minv(a.len() as u32, Tables::MOD.0 as u32) as i64;
    a.iter_mut()
        .for_each(|an| *an = *an * invsize % Tables::MOD.0);
}

fn convolution_mod1(a: &[i64], b: &[i64]) -> Vec<i64> {
    let s = a.len() + b.len() - 1;
    let mut va = NTT_TABLE.0.ntt(a.iter().copied(), s);
    //va.resize(ex, 0);
    let vb = NTT_TABLE.0.ntt(b.iter().copied(), s);
    let ex = va.len();
    //vb.resize(ex, 0);
    va.iter_mut()
        .zip(vb)
        .for_each(|(an, bn)| *an = *an * bn % Tables::MOD.0);
    NTT_TABLE.0.intt_mut_slc_unchecked(va.as_mut_slice());
    let invsize = super::minv(ex as u32, Tables::MOD.0 as u32) as i64;
    va.iter_mut()
        .for_each(|an| *an = *an * invsize % Tables::MOD.0);
    va
}

fn convolution_mod<const M: u32, const G: u32>(
    a: impl Iterator<Item = i64>,
    b: impl Iterator<Item = i64>,
    size_a: usize,
    size_b: usize,
    tbl: &NttTable<M, G>,
) -> Vec<i64> {
    let s = size_a + size_b - 1;
    let mut va = tbl.ntt(a, s);
    let vb = tbl.ntt(b, s);
    let ex = va.len();
    va.iter_mut()
        .zip(vb)
        .for_each(|(an, bn)| *an = *an * bn % M as i64);
    tbl.intt_mut_slc_unchecked(&mut va);
    let invsize = super::minv(ex as u32, M) as i64;
    va.iter_mut().for_each(|an| *an = *an * invsize % M as i64);
    va
}

fn convolution_mod2_mut_unchecked(a: &mut [i64], b: &mut [i64]) {
    NTT_TABLE.1.ntt_mut_slc_unchecked(a);
    NTT_TABLE.1.ntt_mut_slc_unchecked(b);
    a.iter_mut()
        .zip(b.iter())
        .for_each(|(an, bn)| *an = *an * bn % Tables::MOD.1);
    NTT_TABLE.1.intt_mut_slc_unchecked(a);
}

fn convolution_mod2(a: &[i64], b: &[i64]) -> Vec<i64> {
    let s = a.len() + b.len() - 1;
    let ex = s.next_power_of_two();
    let mut va = a
        .iter()
        .map(|an| an.rem_euclid(Tables::MOD.1))
        .chain((0..(ex - a.len())).map(|_| 0))
        .collect::<Vec<i64>>();
    //va.resize(ex, 0);
    let mut vb = b
        .iter()
        .map(|bn| bn.rem_euclid(Tables::MOD.1))
        .chain((0..(ex - b.len())).map(|_| 0))
        .collect::<Vec<i64>>();
    //vb.resize(ex, 0);
    NTT_TABLE.1.ntt_mut_slc_unchecked(va.as_mut_slice());
    NTT_TABLE.1.ntt_mut_slc_unchecked(vb.as_mut_slice());
    va.iter_mut()
        .zip(vb)
        .for_each(|(an, bn)| *an = *an * bn % Tables::MOD.1);
    NTT_TABLE.1.intt_mut_slc_unchecked(va.as_mut_slice());
    let invsize = super::minv(ex as u32, Tables::MOD.1 as u32) as i64;
    va.iter_mut()
        .for_each(|an| *an = *an * invsize % Tables::MOD.1);
    va
}

fn convolution_mod3_mut_unchecked(a: &mut [i64], b: &mut [i64]) {
    NTT_TABLE.2.ntt_mut_slc_unchecked(a);
    NTT_TABLE.2.ntt_mut_slc_unchecked(b);
    a.iter_mut()
        .zip(b.iter())
        .for_each(|(an, bn)| *an = *an * bn % Tables::MOD.2);
    NTT_TABLE.2.intt_mut_slc_unchecked(a);
}

fn convolution_mod3(a: &[i64], b: &[i64]) -> Vec<i64> {
    let s = a.len() + b.len() - 1;
    let ex = s.next_power_of_two();
    let mut va = a
        .iter()
        .map(|an| an.rem_euclid(Tables::MOD.2))
        .chain((0..(ex - a.len())).map(|_| 0))
        .collect::<Vec<i64>>();
    //va.resize(ex, 0);
    let mut vb = b
        .iter()
        .map(|bn| bn.rem_euclid(Tables::MOD.2))
        .chain((0..(ex - b.len())).map(|_| 0))
        .collect::<Vec<i64>>();
    //vb.resize(ex, 0);
    NTT_TABLE.2.ntt_mut_slc_unchecked(va.as_mut_slice());
    NTT_TABLE.2.ntt_mut_slc_unchecked(vb.as_mut_slice());
    va.iter_mut()
        .zip(vb)
        .for_each(|(an, bn)| *an = *an * bn % Tables::MOD.2);
    NTT_TABLE.2.intt_mut_slc_unchecked(va.as_mut_slice());
    let invsize = super::minv(ex as u32, Tables::MOD.2 as u32) as i64;
    va.iter_mut()
        .for_each(|an| *an = *an * invsize % Tables::MOD.2);
    va
}

fn delete_zero(v: &mut Vec<i64>) {
    if let Some(i) = v.iter().rposition(|n| *n != 0) {
        v.truncate(i + 1);
    } else {
        *v = vec![0];
    }
}

/// aとbの畳み込みの結果がすべて2013265921未満でなければならない
fn convolution_level1(
    a: impl Iterator<Item = i64> + Send,
    b: impl Iterator<Item = i64>,
    size_a: usize,
    size_b: usize,
) -> Vec<i64> {
    NTT_TABLE.0.convolution_parallel(a, b, size_a, size_b)
}

pub fn convolution_level2(
    a: impl Iterator<Item = i64> + Send,
    b: impl Iterator<Item = i64>,
    size_a: usize,
    size_b: usize,
) -> Vec<i64> {
    let (mut v1, mut v2) = (&NTT_TABLE.0, &NTT_TABLE.1).convolution_parallel(a, b, size_a, size_b);
    v1.truncate(size_a + size_b - 1);
    v2.truncate(v1.len());
    v1.iter_mut()
        .zip(v2.iter())
        .for_each(|(n1, n2)| *n1 = garner2(*n1, *n2, Tables::MOD.0, Tables::MOD.1));
    v1
}

const fn garner2(a: i64, b: i64, m1: i64, m2: i64) -> i64 {
    a + ((b - a) * super::minv(m1 as u32, m2 as u32) as i64).rem_euclid(m2) * m1
}

#[cfg(test)]
mod test {
    use std::time::Instant;

    #[allow(unused_imports)]
    use super::*;
    use rand::Rng;

    #[test]
    fn convolute() {
        let a = &[2, 2, 20, 2];
        let b = &[1, 1, 1, 1];
        let v = convolution_level2(a.iter().copied(), b.iter().copied(), a.len(), b.len());
        println!("{v:?}");
    }

    #[test]
    fn time() {
        let n = 1 << 22;
        let mut rng = rand::thread_rng();
        let a: Vec<i64> = (0..n).map(|_| rng.gen_range(0..1000000000)).collect();
        let b: Vec<i64> = (0..n).map(|_| rng.gen_range(0..1000000000)).collect();
        let s = Instant::now();
        let v = convolution_level2(a.into_iter(), b.into_iter(), n, n);
        println!("elapsed: {} ms", s.elapsed().as_micros() as f32 / 1000.);
    }
}
