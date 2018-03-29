// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use ndarray::{Array};
use std::f32;
use std::f64;
use statrs::distribution::{ChiSquared, Univariate, Continuous};

fn log_sum_exp(p: &Vec<f64>) -> f64{
    let max_p: f64 = p.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let sum_rst: f64 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln()
}

pub struct SmoothingInfo {
    pub num_states:     usize,
    pub num_positions:  usize,
    pub forward_prob:   Vec<Vec<f64>>,
    pub backward_prob:  Vec<Vec<f64>>,
    pub status_prob:    Vec<Vec<f64>>, // prob(x_n | o_1, o_2, ......, o_N)
    pub null_state:     usize,
    pub compared_to_second_best:    bool,
    pub max_likely_state:   Vec<u16>,
    pub max_loglikely:  Vec<f64>,
    pub pval:           Vec<f64>,
    pub max_likely_sequence:    Vec<u16>,
}

pub trait HmmModel {
    fn num_positions(&self) -> usize;
    fn num_states(&self) -> usize;

    fn prior(&self, state: usize) -> f64;
    fn score_state(&self, pos: usize, state: usize) -> f64;
    fn score_trans(&self, pos: usize, prev_state: usize, state: usize) -> f64;

    fn viterbi(&self) -> Vec<u16> {
        let mut tb = Array::from_elem((self.num_positions(), self.num_states()), 0);
        let mut vit = Array::from_elem((self.num_positions(), self.num_states()), f32::NEG_INFINITY);

        // Copy over prior
        for i in 0..self.num_states() {
            vit[(0,i)] = self.prior(i) as f32;
        }

        for i in 1..self.num_positions() {

            for state in 0..self.num_states() {
                let state_score = self.score_state(i, state) as f32;

                for prev_state in 0..self.num_states() {

                    let trans_score = self.score_trans(i, prev_state, state) as f32;
                    let score = vit[(i-1, prev_state)] + trans_score + state_score;

                    if score > vit[(i, state)]
                    {
                        vit[(i,state)] = score;
                        tb[(i,state)] = prev_state;
                    }
                }
            }
        }

        // Find best final state
        let mut best_final_state = 0;
        let mut best_final_score = f32::NEG_INFINITY;
        let last_pos = self.num_positions() - 1;

        for st in 0..self.num_states() {
            if vit[(last_pos, st)] > best_final_score {
                best_final_state = st;
                best_final_score = vit[(last_pos, st)];
            }
        }

        let mut state_vec = Array::from_elem((self.num_positions()), 0 as u16);
        let mut pos = last_pos;
        let mut current_state = best_final_state;

        loop {
            state_vec[pos] = current_state as u16;
            current_state = tb[(pos, current_state)];

            if pos == 0 {
                break;
            }
            pos = pos - 1;
        }

        let cp: Vec<u16> = state_vec.iter().cloned().collect();
        cp
    }

    fn forward(&self) -> Vec<Vec<f64>> {
        let ns = self.num_states();
        let np = self.num_positions();
        let mut prob = vec![vec![0.0f64; ns]; np];
        for s in 0..ns {
            prob[0][s] = self.prior(s) + self.score_state(0, s);
        }

        for p in 1..np {
            for s in 0..ns {
                prob[p][s] = self.score_state(p, s) + log_sum_exp(&
                    (0..ns).map(|x| prob[p-1][x]+self.score_trans(p, x, s))
                    .collect::<Vec<f64>>());
            }
        }

        prob
    }

    fn backward(&self) -> Vec<Vec<f64>> {
        let ns = self.num_states();
        let np = self.num_positions();
        let mut prob = vec![vec![0.0f64; ns]; np];
        for s in 0..ns {
            prob[np-1][s] = 0.0;
        }
        
        for p in (1..np).rev() {
            for s in 0..ns {
                prob[p-1][s] = log_sum_exp(& (0..ns).map(|x| 
                    prob[p][x] + self.score_state(p, x) + self.score_trans(p, s, x))
                    .collect::<Vec<f64>>());
            }
        }
        prob
    }

    fn get_inference(&self, default_state: u16) -> SmoothingInfo {
        let ns = self.num_states();
        let np = self.num_positions();
        SmoothingInfo {
            num_states:     ns,
            num_positions:  np,
            forward_prob:   self.forward(),
            backward_prob:  self.backward(),
            status_prob:    vec![vec![0.0f64; ns]; np],
            null_state:     default_state as usize,
            max_likely_state:   vec![0u16; np],
            max_loglikely:      vec![0.0f64; np],
            pval:           vec![0.0f64; np],
            max_likely_sequence:    self.viterbi(),
            compared_to_second_best: false,
        }
    }
}

impl SmoothingInfo {
    pub fn compute(&mut self) {
        // calculate the status probability prob(x_n | o_1, o_2, ......, o_N) 
        let chi_sq = ChiSquared::new(1.0).unwrap();
        for p in 0..self.num_positions {
            let mut max_ll = f64::NEG_INFINITY;
            let mut sec_ll = f64::NEG_INFINITY;
            let mut best_state = 0usize;
            for s in 0..self.num_states {
                self.status_prob[p][s] = self.forward_prob[p][s] + self.backward_prob[p][s];
                if self.status_prob[p][s] > max_ll {
                    sec_ll = max_ll;
                    max_ll = self.status_prob[p][s];
                    best_state = s;
                } else if self.status_prob[p][s] > sec_ll {
                    sec_ll = self.status_prob[p][s];
                }
            }

            self.max_likely_state[p] = best_state as u16;
            self.max_loglikely[p] = (1.0e-60f64).max(2.0 * (max_ll  - if self.compared_to_second_best { sec_ll } else {self.status_prob[p][self.null_state]})); 
            self.pval[p] = 1.0 - chi_sq.cdf(self.max_loglikely[p])+chi_sq.pdf(self.max_loglikely[p]);
        }
    }

    pub fn set_compared_to_second_best(&mut self, use_second_best: bool) {
        self.compared_to_second_best = use_second_best;
    }

    pub fn interval_significance(&self, hmm: &HmmModel, start_pos: usize, end_pos: /*non inclusive*/ usize) -> (u16 /*state*/,f64 /*pval*/, f64 /*lld*/) {
        let chi_sq = ChiSquared::new(1.0).unwrap();
        let mut max_ll = f64::NEG_INFINITY;
        let mut sec_ll = f64::NEG_INFINITY;
        let mut best_state = 0usize;
        let mut ll = vec![0.0f64; self.num_states];
        for s in 0..self.num_states {
            ll[s] = self.forward_prob[start_pos][s] + self.backward_prob[end_pos-1][s];
            for p in start_pos+1..end_pos {
                ll[s] += hmm.score_state(p, s) + hmm.score_trans(p, s, s);
            }
            if ll[s] > max_ll {
                sec_ll = max_ll;
                max_ll = ll[s];
                best_state = s;
            } else if ll[s] > sec_ll {
                sec_ll = ll[s];
            }
        }


        let lld = (1.0e-60f64).max(2.0 * (max_ll  - if self.compared_to_second_best { sec_ll } else {ll[self.null_state]})); 
        let mut pval = 1.0-chi_sq.cdf(lld)+chi_sq.pdf(lld);
        if pval > 1.0 {
            pval = 1.0;
        } else if pval < 0.0 {
            pval = 0.0;
        }

        println!("best_state: {} start: {} end: {}\t lld {:.8e}\tpval {:.8e}", best_state, start_pos, end_pos, lld, pval);
        (best_state as u16, lld, pval)
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use statrs::distribution::{ChiSquared, Univariate};

    struct ObsVec {
        p_error: f64,
        obs: Vec<u16>
    }

    impl HmmModel for ObsVec {
        fn num_positions(&self) -> usize {
            self.obs.len()
        }

        fn num_states(&self) -> usize {
            2
        }

        fn prior(&self, _: usize) -> f64 {
            0.5f64.ln()
        }

        fn score_state(&self, pos: usize, state: usize) -> f64 {
            let obs = self.obs[pos];
            if obs == (state as u16) {
                (1.0 - self.p_error).ln()
            } else {
                self.p_error.ln()
            }
        }

        #[allow(unused_variables)]
        fn score_trans(&self, pos: usize, prev_state: usize, state: usize) -> f64 {
            if prev_state == state {
                0.8f64.ln()
            } else {
                0.2f64.ln()
            }
        }
    }

    fn simple_example() -> ObsVec {
        let v : Vec<u16> = vec![1,1,1,0,1,1,1,0,0,0,0,1,0,0,0,0];
        ObsVec {
            p_error: 0.25,
            obs: v
        }
    }

    #[test]
    fn simple_switch() {
        let obs = simple_example();

        let traceback = obs.viterbi();
        assert_eq!(traceback, vec![1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0]);
    }

    #[test]
    fn test_forwad() {
        let obs = simple_example();
        let fwd_prob = obs.forward();
        let exp = vec![
            vec![1.2500e-1,       3.7500e-1],
            vec![4.3750e-2,       2.4375e-1],
            vec![2.0938e-2,       1.5281e-1],
            vec![3.5484e-2,       3.1609e-2],
            vec![8.6773e-3,       2.4288e-2],
            vec![2.9499e-3,       1.5875e-2],
            vec![1.3837e-3,       9.9672e-3],
            vec![2.3253e-3,       2.0626e-3],
            vec![1.7046e-3,       5.2879e-4],
            vec![1.1021e-3,       1.9099e-4],
            vec![6.8989e-4,       9.3301e-5],
            vec![1.4264e-4,       1.5946e-4],
            vec![1.0951e-4,       3.9025e-5],
            vec![7.1557e-5,       1.3280e-5],
            vec![4.4926e-5,       6.2339e-6],
            vec![2.7891e-5,       3.4931e-6],
        ];
        for i in 0..obs.num_positions() {
            let r = fwd_prob[i].iter().map(|x| x.exp()).collect::<Vec<f64>>();
            assert!( (r[0]-exp[i][0]).abs() < 1e-4 && (r[1]-exp[i][1]).abs() < 1e-4, "not the same\nobs {:?}\nexp {:?}\n", r, exp[i]);
        }
    }

    #[test]
    fn test_backwad() {
        let obs = simple_example();
        println!("number of positions {}", obs.num_positions());

        let bwd_prob = obs.backward();
        let exp = vec![
            vec![3.08730732e-5,        7.33991628e-5],
            vec![6.67908399e-5,        1.16766035e-4],
            vec![2.00529767e-4,        1.77899244e-4],
            vec![2.77431032e-4,        6.81422946e-4],
            vec![5.71068244e-4,        1.08811589e-3],
            vec![1.59487612e-3,        1.68062014e-3],
            vec![6.26517911e-3,        2.27893531e-3],
            vec![1.01252361e-2,        3.80074949e-3],
            vec![1.63111977e-2,        6.77034922e-3],
            vec![2.59886406e-2,        1.43602656e-2],
            vec![3.98196875e-2,        4.19365625e-2],
            vec![1.56456250e-1,        5.68562500e-2],
            vec![2.52875000e-1,        9.46250000e-2],
            vec![4.07500000e-1,        1.67500000e-1],
            vec![6.50000000e-1,        3.50000000e-1],
            vec![1.00000000e00,        1.00000000e00],
        ];
        for i in 0..obs.num_positions() {
            let r = bwd_prob[i].iter().map(|x| x.exp()).collect::<Vec<f64>>();
            println!("{}: {:.8e}\t{:.8e}", i, r[0], r[1]);
            assert!( (r[0]-exp[i][0]).abs() < 1e-6 && (r[1]-exp[i][1]).abs() < 1e-6, "not the same\nobs {:?}\nexp {:?}\n", r, exp[i]);
        }
    }

    #[test]
    fn test_inference() {
        let obs = simple_example();
        let mut inference = obs.get_inference(0u16);
        inference.set_compared_to_second_best(true);
        inference.compute();
        for p in 0..inference.num_positions {
            println!("Position {}\tbest state {}\tlld {:.4e}\tp-value {:.4e}\tmax seq {}",
                     p, inference.max_likely_state[p], inference.max_loglikely[p], inference.pval[p], inference.max_likely_sequence[p]);
        }

        println!("\n\n");
        println!("{:?}", inference.interval_significance(&obs, 0,1));
        println!("{:?}", inference.interval_significance(&obs, 0,7));
        println!("{:?}", inference.interval_significance(&obs, 8,16));
        println!("{:?}", inference.interval_significance(&obs, 5,9));

    }
}
    

