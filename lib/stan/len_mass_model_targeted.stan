functions {

    // Interpolate a value at a postion in vals, which has a fixed bin_size
    real interp(int pos, int bin_size, int[] vals, int L) {
        int b1;
        int b2;
        real f;

        b1 <- min(L-1, pos/bin_size + 1);
        b2 <- b1 + 1;

        f <- ((1.0) * pos - ((pos / bin_size)*bin_size)) / bin_size;
        return vals[b2] * f + vals[b1] * (1.0 - f);
    }

    // Approximate the number of target bases in the interval
    // [pos, pos+len) use the discretized cumulative target bases track
    real target_size(int pos, int len, int bin_size, int[] cum_target_bases, int L) {
        real c1;
        real c2;

        c1 <- interp(pos, bin_size, cum_target_bases, L);
        c2 <- interp(pos + len, bin_size, cum_target_bases, L);

        return c2 - c1;
    }
}


data {
    int<lower=1> BC;  // number of BCS
    int<lower=1> N;   // number of fragments
    int<lower=10> K;  // number of length bins

    int<lower=1> bc_observed_frags[BC];        // Observed fragments per BC
    int<lower=1> bin_length[K];               // Fragment length bins
    int<lower=1, upper=BC> bc_id[N];           // BC id of each fragment
    int<lower=0> num_reads[N];                 // Num reads on each fragment
    int<lower=0> obs_length[N];               // Observed length of fragment
    int<lower=0> pos[N];                      // Left hand position of fragment, decimated 8x

    int<lower=1> GB;                           // Genome bins
    int<lower=1> gb_size;                      // Genome bin size in bases, decimated 8x
    int<lower=0> cum_target_bases[GB];         // cum_target_bases[x] is sum_{i = 0 -> gb_size*x} OnTarget(i)
    real<lower=1.0> genome_size;                 // full genome size
}

transformed data {
    vector[K] theta_prior;

    // Uninformative prior for Dirichlet distribution
    theta_prior <- rep_vector(1.05, K);
}


parameters {
    vector<lower=0.0>[BC] alpha;                  // BC amp rate (reads / base pair)
    vector<lower=0>[BC] bc_unobserved_frags;    // # of unobserved fragments per BC
    simplex[K] theta;                           // Fragment length distribution
    real<lower=0> mean_frags;                   // Mean number of fragments in a BC
    real<lower=1> read_disp;
}


transformed parameters {
    real fl_lik[N,K];
    real unobs_lik[BC,K];
    real bc_total_frags[BC];
    real log_theta[K];
    vector[K-1] diff;
    vector[K-2] diff2;

    diff <- head(theta, K-1) - tail(theta, K-1);
    diff2 <- head(diff, K-2) - tail(diff, K-2);

    for (k in 1:K)
        log_theta[k] <- log(theta[k]);

    // Likelihood of the observed num reads and length of each read, given each fragment length
    for (n in 1:N) 
    {
        for (k in 1:K)
        {
            if (bin_length[k] < obs_length[n])
            {
                fl_lik[n,k] <- -100000;
            }
            else 
            {
                int slack;
                int slack_dec;
                int init_win;
                int slack_win;
                int slack_win_dec;

                int slack_bins;
                real target_bases;

                # Determine possible fragment start points
                slack <- (bin_length[k] - obs_length[n]);
                if (slack/8 > pos[n])
                {
                    slack <- pos[n] * 8;
                }

                slack_bins <- min(16, slack / 2048 + 1);
                slack_win <- slack / slack_bins;
                slack_win_dec <- slack / (slack_bins * 8);
                slack_dec <- slack / 8;

                {
                    vector[slack_bins] slack_sum;
                    # Integrate over the slack
                    for (s in 1:slack_bins) {
                        target_bases <- bin_length[k] * 0.01 + fmax(0.0, target_size(pos[n] - slack_dec + slack_win_dec/2 + slack_win_dec*(s-1), bin_length[k]/8, gb_size/8, cum_target_bases, GB));
                        
                        # TODO -- use simplified version of this expr. for speed 
                        #slack_sum[s] <- poisson_log(num_reads[n], alpha[bc_id[n]] * target_bases) - num_reads[n] * log(target_bases);
                        slack_sum[s] <- neg_binomial_2_log(num_reads[n], alpha[bc_id[n]] * target_bases, read_disp) - num_reads[n] * log(target_bases);
                    }

                    fl_lik[n,k] <- log_theta[k] + log((1.0 * slack_win) / genome_size) + log_sum_exp(slack_sum);
                }
            }
        }
    }

    // Total number of fragments present on BC
    for (b in 1:BC)
        bc_total_frags[b] <- bc_observed_frags[b] + bc_unobserved_frags[b];

    // Likelihood of an unobserved fragment in BC b, given a fragment length
    for (b in 1:BC)
        for (k in 1:K)
            #unobs_lik[b,k] <- log(theta[k]) + poisson_log(0, bin_length[k] * alpha[b]);
            unobs_lik[b,k] <- log(theta[k]) + neg_binomial_2_log(0, bin_length[k] * alpha[b] * 0.1, read_disp);
}


model {
    // Fragment length distribution
    theta ~ dirichlet(theta_prior);
    increment_log_prob(normal_log(diff2, 0.0, 0.5/K));

    // Likelihood of the number of fragments in a BC -- use continuous analog of Poisson
    increment_log_prob(gamma_log(bc_total_frags, mean_frags, 1.0));

    // Likelihood of unobserved fragments
    for (b in 1:BC)
        increment_log_prob(bc_unobserved_frags[b] * log_sum_exp(unobs_lik[b]) + binomial_coefficient_log(bc_total_frags[b], bc_unobserved_frags[b]));

    // Likelihood of observed fragments: sum over fragment length mixture, product over fragments
    for (n in 1:N)
        increment_log_prob(log_sum_exp(fl_lik[n]));
}
