pub mod convolution;
use std::{iter::Iterator, num::NonZeroUsize};

#[derive(Debug)]
struct NttTable<const M: u32, const G: u32> {
    root1: [u32; 32],
    iroot1: [u32; 32],
    root3: [u32; 30],
    iroot3: [u32; 30],
}

const fn powm<const M: u32>(base: u32, exp: u32) -> u32 {
    let mut mask = exp.next_power_of_two();
    let mut ans = 1;
    while mask > 0 {
        ans = (ans as u64 * ans as u64 % M as u64) as u32;
        if exp & mask != 0 {
            ans = (ans as u64 * base as u64 % M as u64) as u32;
        }
        mask >>= 1;
    }
    ans
}

const fn minv(x: u32, m: u32) -> u32 {
    let (mut a, mut b) = (x as i64, m as i64);
    let (mut u, mut v) = (1, 0);
    while b != 0 {
        let tmp = u;
        u = v;
        v = tmp - (a / b) * v;

        let tmp = a;
        a = b;
        b = tmp % b;
    }
    u %= m as i64;
    if u.is_negative() {
        u += m as i64;
    }
    u as u32
}

impl<const M: u32, const G: u32> NttTable<M, G> {
    const fn new() -> Self {
        let rank = (M - 1).trailing_zeros() as usize;
        let mut root1 = [0; 32];
        let mut iroot1 = [0; 32];
        root1[rank] = powm::<M>(G, M >> rank as u32);
        iroot1[rank] = minv(root1[rank], M);
        let mut i = rank;
        while i > 0 {
            i -= 1;
            root1[i] = (root1[i + 1] as u64 * root1[i + 1] as u64 % M as u64) as u32;
            iroot1[i] = (iroot1[i + 1] as u64 * iroot1[i + 1] as u64 % M as u64) as u32;
        }
        let mut root3 = [0; 30];
        let mut iroot3 = [0; 30];
        i = rank - 2;
        while i > 0 {
            i -= 1;
            root3[i] = (iroot1[2] as u64 * root1[i + 2] as u64 % M as u64 * root1[i + 3] as u64
                % M as u64) as u32;
            iroot3[i] = (root1[2] as u64 * iroot1[i + 2] as u64 % M as u64 * iroot1[i + 3] as u64
                % M as u64) as u32;
        }
        Self {
            root1,
            iroot1,
            root3,
            iroot3,
        }
    }

    // sの長さが2のべき乗であることを仮定
    fn ntt_mut_slc_unchecked(&self, s: &mut [i64]) {
        let n = s.len();
        let mut order = n.trailing_zeros();
        let m2 = M as u64 * M as u64;
        let w_4 = self.root1[2];
        while order > 1 {
            let mut w = 1;
            let loops = 1 << (order - 2);
            for i in 0..(n >> order) {
                let offset = i << order;
                let w2 = (w as u64 * w as u64 % M as u64) as u32;
                let w3 = (w as u64 * w2 as u64 % M as u64) as u32;
                for j in 0..loops {
                    let k = j + offset;
                    let a1 = s[k] as u64;
                    let a2 = w as u64 * s[k + loops] as u64;
                    let a3 = w2 as u64 * s[k + loops * 2] as u64;
                    let a4 = w3 as u64 * s[k + loops * 3] as u64;
                    let a2a4 = w_4 as u64 * ((a2 + m2 - a4) % M as u64);
                    s[k] = ((a1 + a2 + a3 + a4) % M as u64) as i64;
                    s[k + loops] = ((a1 + a3 + m2 * 2 - (a2 + a4)) % M as u64) as i64;
                    s[k + loops * 2] = ((a1 + m2 - a3 + a2a4) % M as u64) as i64;
                    s[k + loops * 3] = ((a1 + m2 * 2 - a3 - a2a4) % M as u64) as i64;
                }
                w = (w as u64 * self.root3[i.trailing_ones() as usize] as u64 % M as u64) as u32;
            }
            order -= 2;
        }
        if order == 1 {
            let mut w = 1;
            for i in 0..(n / 4) {
                let j = i * 4;
                let a1 = s[j] as u64;
                let a2 = w as u64 * s[j + 1] as u64;
                let a3 = s[j + 2] as u64;
                let a4 = w as u64 * s[j + 3] as u64 % M as u64 * w_4 as u64;
                s[j] = ((a1 + a2) % M as u64) as i64;
                s[j + 1] = ((a1 + m2 - a2) % M as u64) as i64;
                s[j + 2] = ((a3 + a4) % M as u64) as i64;
                s[j + 3] = ((a3 + m2 - a4) % M as u64) as i64;
                w = (w as u64 * self.root3[i.trailing_ones() as usize] as u64 % M as u64) as u32;
            }
        }
    }

    fn intt_mut_slc_unchecked(&self, s: &mut [i64]) {
        let n = s.len();
        let o = n.trailing_zeros();
        let mut order = o;
        let w_4 = self.iroot1[2];
        if (order & 1) != 0 {
            let mut w = 1;
            for i in 0..(n / 4) {
                let j = i * 4;
                let a1 = s[j] as u32;
                let a2 = s[j + 1] as u32;
                let a3 = s[j + 2] as u32;
                let a4 = s[j + 3] as u32;
                s[j] = ((a1 + a2) % M) as i64;
                s[j + 1] = (w as u64 * (a1 + M - a2) as u64 % M as u64) as i64;
                s[j + 2] = ((a3 + a4) % M) as i64;
                s[j + 3] =
                    (w as u64 * w_4 as u64 % M as u64 * (a3 + M - a4) as u64 % M as u64) as i64;
                w = (w as u64 * self.iroot3[i.trailing_ones() as usize] as u64 % M as u64) as u32;
            }
            order -= 1;
        }
        while order > 0 {
            let mut w = 1;
            let loops = n >> order;
            for i in 0..(1 << (order - 2)) {
                let offset = i << (o - order + 2);
                let w2 = (w as u64 * w as u64 % M as u64) as u32;
                let w3 = (w as u64 * w2 as u64 % M as u64) as u32;
                for j in 0..loops {
                    let k = j + offset;
                    let a1 = s[k] as u64;
                    let a2 = s[k + loops] as u64;
                    let a3 = s[k + loops * 2] as u64;
                    let a4 = s[k + loops * 3] as u64;
                    let a3a4 = w_4 as u64 * (a3 + M as u64 - a4) % M as u64;
                    s[k] = ((a1 + a2 + a3 + a4) % M as u64) as i64;
                    s[k + loops] = (w as u64 * (a1 + M as u64 - a2 + a3a4) % M as u64) as i64;
                    s[k + loops * 2] =
                        (w2 as u64 * (a1 + a2 + M as u64 * 2 - a3 - a4) % M as u64) as i64;
                    s[k + loops * 3] =
                        (w3 as u64 * (a1 + M as u64 * 2 - a2 - a3a4) % M as u64) as i64;
                }
                w = (w as u64 * self.root3[i.trailing_ones() as usize] as u64 % M as u64) as u32;
            }
            order -= 2;
        }
    }

    fn ntt(&self, i: impl Iterator<Item = i64>, size: usize) -> Vec<i64> {
        let mut v: Vec<i64> = i.take(size).map(|n| n.rem_euclid(M as i64)).collect();
        v.resize(size.next_power_of_two(), 0);
        self.ntt_mut_slc_unchecked(&mut v);
        v
    }

    fn intt(&self, i: impl Iterator<Item = i64>, size: usize) -> Vec<i64> {
        let mut v: Vec<i64> = i.take(size).map(|n| n.rem_euclid(M as i64)).collect();
        v.resize(size.next_power_of_two(), 0);
        self.intt_mut_slc_unchecked(&mut v);
        v
    }
}

pub trait ParallelConvoluter {
    type Ret;
    fn convolution_parallel(
        self,
        a: impl Iterator<Item = i64> + Send,
        b: impl Iterator<Item = i64>,
        size_a: usize,
        size_b: usize,
    ) -> Self::Ret;
}

impl<'a, const M: u32, const G: u32> ParallelConvoluter for &'a NttTable<M, G> {
    type Ret = Vec<i64>;
    fn convolution_parallel(
        self,
        a: impl Iterator<Item = i64> + Send,
        b: impl Iterator<Item = i64>,
        size_a: usize,
        size_b: usize,
    ) -> Self::Ret {
        let size = size_a + size_b - 1;
        let (vb, mut va) = thread::scope(|s| {
            let handle = s.spawn(|| self.ntt(a, size));
            (self.ntt(b, size), handle.join().unwrap())
        });
        let ex = va.len();
        va.iter_mut()
            .zip(vb)
            .for_each(|(an, bn)| *an = *an * bn % M as i64);
        self.intt_mut_slc_unchecked(&mut va);
        let invsize = minv(ex as u32, M) as i64;
        va.iter_mut().for_each(|an| *an = *an * invsize % M as i64);
        va
    }
}
use std::thread;

impl<'a, 'b, const M1: u32, const G1: u32, const M2: u32, const G2: u32> ParallelConvoluter
    for (&'a NttTable<M1, G1>, &'b NttTable<M2, G2>)
{
    type Ret = (Vec<i64>, Vec<i64>);
    fn convolution_parallel(
        self,
        a: impl Iterator<Item = i64> + Send,
        b: impl Iterator<Item = i64>,
        size_a: usize,
        size_b: usize,
    ) -> Self::Ret {
        let s = size_a + size_b - 1;
        let mut va: (Vec<i64>, Vec<i64>) = a
            .take(s)
            .map(|n| (n.rem_euclid(M1 as i64), n.rem_euclid(M2 as i64)))
            .unzip();
        let mut vb: (Vec<i64>, Vec<i64>) = b
            .take(s)
            .map(|n| (n.rem_euclid(M1 as i64), n.rem_euclid(M2 as i64)))
            .unzip();
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
                    .for_each(|(an, bn)| *an = *an * bn % M2 as i64);
                let invsize = minv(ex as u32, M2) as i64;
                va.1.iter_mut()
                    .for_each(|an| *an = *an * invsize % M2 as i64);
                self.1.intt_mut_slc_unchecked(&mut va.1);
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
                .for_each(|(an, bn)| *an = *an * bn % M1 as i64);
            let invsize = minv(ex as u32, M1) as i64;
            va.0.iter_mut()
                .for_each(|an| *an = *an * invsize % M1 as i64);
            self.0.intt_mut_slc_unchecked(&mut va.0);
        });
        (va.0, va.1)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::Rng;
    use std::time::Instant;

    #[test]
    fn nttnew() {
        let tbl = NttTable::<17, 3>::new();
        println!("{:?}", tbl);
    }

    #[test]
    fn nttmutslcunchecked() {
        let tbl = NttTable::<17, 3>::new();
        let slc = &mut [1, 2, 3, 4];
        let slc2 = &mut [1, 1, 1, 1];
        tbl.ntt_mut_slc_unchecked(slc);
        tbl.ntt_mut_slc_unchecked(slc2);
        println!("slc: {:?}\n slc2:{:?}", slc, slc2);
    }

    #[test]
    fn overflowcheck() {
        let tbl = NttTable::<2013265921, 31>::new();
        let mut rng = rand::thread_rng();
        let n = 1 << 20;
        let mut v = Vec::with_capacity(n);
        v.resize_with(n, || rng.gen_range(0..2013265921));
        let s = Instant::now();
        tbl.ntt_mut_slc_unchecked(v.as_mut_slice());
        println!("elapsed: {} ms", s.elapsed().as_micros() as f32 / 1000.0f32);
    }

    #[test]
    fn inversecheck() {
        let tbl = NttTable::<1811939329, 13>::new();
        let slc = &mut [3, 4, 5, 6, 1, 3, 5, 7];
        tbl.ntt_mut_slc_unchecked(slc);
        tbl.intt_mut_slc_unchecked(slc);
        println!("{:?}", slc);
    }
}
