data {
    int<lower=1> BC;  // number of BCS
    int<lower=1> N;   // number of fragments
    int<lower=10> K;  // number of length bins

    int<lower=1> bc_observed_frags[BC];        // Observed fragments per BC
    real<lower=1> bin_length[K];               // Fragment length bins
    int<lower=1, upper=BC> bc_id[N];           // BC id of each fragment
    int<lower=0> num_reads[N];                 // Num reads on each fragment
    real<lower=0> obs_length[N];               // Observed length of fragment 
}

transformed data {
    vector[K] theta_prior;
    real k2[N,K];
    real f;
    real frag_l;
    int nr;
    real obs_l;
    real len_diff;

    // Compute the likelihood of the observed length of each fragment, given each fragment length
    for (n in 1:N) {
        nr <- num_reads[n];
        obs_l <- obs_length[n];

        for (k in 1:K) {
            frag_l <- bin_length[k];
            f <- obs_l / frag_l;
            k2[n,k] <- if_else(nr < 2, 0.0, 
                           if_else(bin_length[k] <= obs_length[n], -100000, 
                               log(nr * (nr-1) * (1.0 - f) * (1.0/frag_l)) + (nr-2) * log(f)));
        }
    }

    // Uninformative prior for Dirichlet distribution
    theta_prior <- rep_vector(1.05, K);
}

parameters {
    vector<lower=0.0>[BC] alpha;                  // BC amp rate (reads / base pair)
    vector<lower=0>[BC] bc_unobserved_frags;    // # of unobserved fragments per BC
    simplex[K] theta;                           // Fragment length distribution
    real<lower=1> mean_frags;                   // Mean number of fragments in a BC
    real<lower=1> read_disp;
    real<lower=0> amp_length_k;
}

transformed parameters {
    real fl_lik[N,K];
    real unobs_lik[BC,K];
    real bc_total_frags[BC];

    vector[K-1] diff;
    vector[K-2] diff2;
    diff <- head(theta, K-1) - tail(theta, K-1);
    diff2 <- head(diff, K-2) - tail(diff, K-2);

    // Likelihood of the observed num reads and length of each read, given each fragment length
    for (n in 1:N) 
        for (k in 1:K) 
            fl_lik[n,k] <- log(theta[k]) + 
                           neg_binomial_2_log(num_reads[n], (bin_length[k] + bin_length[k]^2 * amp_length_k) * alpha[bc_id[n]], read_disp) + k2[n,k];

    // Total number of fragments present on BC
    for (b in 1:BC)
        bc_total_frags[b] <- bc_observed_frags[b] + bc_unobserved_frags[b];

    // Likelihood of an unobserved fragment in BC b, given a fragment length
    for (b in 1:BC)
        for (k in 1:K)
            unobs_lik[b,k] <- log(theta[k]) + neg_binomial_2_log(0, (bin_length[k] + bin_length[k]^2 * amp_length_k) * alpha[b], read_disp);
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
