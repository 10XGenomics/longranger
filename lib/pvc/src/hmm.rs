// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use ndarray::Array;
use std::f32;

pub trait HmmModel {
    fn num_positions(&self) -> usize;
    fn num_states(&self) -> usize;

    fn prior(&self, state: usize) -> f64;
    fn score_state(&self, pos: usize, state: usize) -> f64;
    fn score_trans(&self, pos: usize, prev_state: usize, state: usize) -> f64;
}

pub fn viterbi<M: HmmModel>(model: M) -> Vec<u16> {
    let mut tb = Array::from_elem((model.num_positions(), model.num_states()), 0);
    let mut vit = Array::from_elem((model.num_positions(), model.num_states()),
                                   f32::NEG_INFINITY);

    // Copy over prior
    for i in 0..model.num_states() {
        vit[(0, i)] = model.prior(i).ln() as f32;
    }

    for i in 1..model.num_positions() {

        for state in 0..model.num_states() {
            let state_score = model.score_state(i, state).ln() as f32;

            for prev_state in 0..model.num_states() {

                let trans_score = model.score_trans(i, prev_state, state).ln() as f32;
                let score = vit[(i - 1, prev_state)] + trans_score + state_score;

                if score > vit[(i, state)] {
                    vit[(i, state)] = score;
                    tb[(i, state)] = prev_state;
                }
            }
        }
    }

    //println!("\n{:?}", vit);
    // Find best final state
    let mut best_final_state = 0;
    let mut best_final_score = f32::NEG_INFINITY;
    let last_pos = model.num_positions() - 1;

    for st in 0..model.num_states() {
        if vit[(last_pos, st)] > best_final_score {
            best_final_state = st;
            best_final_score = vit[(last_pos, st)];
        }
    }

    let mut state_vec = Array::from_elem((model.num_positions()), 0 as u16);
    let mut pos = last_pos;
    let mut current_state = best_final_state;

    loop {
        state_vec[pos] = current_state as u16;
        current_state = tb[(pos, current_state)];

        if pos == 0 {
            break;
        }
        pos -= 1;
    }

    state_vec.iter().cloned().collect()
}


#[cfg(test)]
mod tests {
    use super::*;

    struct ObsVec {
        p_error: f64,
        obs: Vec<u16>,
    }

    impl HmmModel for ObsVec {
        fn num_positions(&self) -> usize {
            self.obs.len()
        }

        fn num_states(&self) -> usize {
            2
        }

        fn prior(&self, _: usize) -> f64 {
            0.5
        }

        fn score_state(&self, pos: usize, state: usize) -> f64 {
            let obs = self.obs[pos];
            if obs == (state as u16) {
                1.0 - self.p_error
            } else {
                self.p_error
            }
        }

        fn score_trans(&self, _: usize, prev_state: usize, state: usize) -> f64 {
            if prev_state == state { 0.8 } else { 0.2 }
        }
    }

    #[test]
    fn simple_switch() {

        let v: Vec<u16> = vec![1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0];
        let model = ObsVec {
            p_error: 0.25,
            obs: v,
        };

        let traceback = viterbi(model);
        assert_eq!(traceback,
                   vec![1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
    }
}
