#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Reports length and mass of input DNA
#
import numpy as np
import scipy.stats
import tenkit.hdf5
import tenkit.stats as tk_stats
import tenkit.safe_json
import tenkit.pandas as p
import os
import subprocess
import math
import shutil
import tenkit.bio_io as tk_io
import tenkit.seq as tk_seq
from tenkit.constants import FRAG_LEN_HIST_BIN_SIZE, FLUIDICS_PARAMS
import tenkit.reference
import martian

__MRO__ = '''
call REPORT_LENGTH_MASS(
    in h5  barcodes,
    in h5 fragments,
    in string reference_path,
    in bed targets_file,
    in string barcode_whitelist,
    out json inferred_length_distribution,
    out json summary,
) split using (
)
'''

MIN_BCS = 8
NUM_FRAGS = 30000
NUM_FRAGS_TARGETED = 15000
NUM_LENGTH_BINS = 65

MIN_RPF = 1.4
MIN_TARGET_SIZE = 2e6


def fmt_val_stan(name, val):
    ''' Format a named value (scalar or array) to stan format '''
    if type(val) is list or hasattr(val, '__iter__'):
        return "%s <- c(%s)" % (name, ", ".join([str(v) for v in val]))
    else:
        return "%s <- %s" % (name, repr(val))

def write_stan_input(fn, vdict):
    ''' Write a dictionary of named values to a R 'dump file' for consumption by stan '''
    with open(fn, 'w') as f:
        for name, val in vdict.iteritems():
            f.write(fmt_val_stan(name, val) + "\n")

def gen_bin_length(n_bins, min_len = 500, max_len = 175000):
    ''' Generate the length bins used for fragment length inference '''
    bin_length = np.linspace(min_len, max_len, n_bins).astype(np.int32)
    return {'bin_length': list(bin_length), 'K': len(bin_length)}

def load_stan_output(fn):
    f = open(fn)

    cols = None
    vals = None

    for l in f:
        if l[0] == '#':
            continue
        v = l.strip().split(',')
        if cols is None:
            cols = v
            continue
        if vals is None:
            vals = [float(i) for i in v]
            continue

    vs = {}
    for (k,v) in zip(cols, vals):
        if k.find(".") > 0:
            k_parts = k.split(".")
            coords = tuple(int(x) for x in k_parts[1:])
            arr_val = vs.setdefault(k_parts[0], {})
            arr_val[coords] = v
        else:
            vs[k] = v

    def proc_item(i):
        if type(i) == dict:
            return mk_array(i)

        return i

    return {k : proc_item(v) for (k,v) in vs.items()}

def mk_array(d):
    dims = len(d.keys()[0])
    shape = [0]*dims
    for dim in range(dims):
        shape[dim] = max(k[dim] for k in d.keys())
    arr = np.zeros(tuple(shape))
    for (k,v) in d.iteritems():
        arr[tuple(x-1 for x in k)] = v

    return arr


def split(args):
    mem_gb = 24
    chunk_defs = [{'__mem_gb': mem_gb}]
    return {'chunks': chunk_defs}


def join(args, outs, chunk_defs, chunk_outs):
    chunk_out = chunk_outs[0]
    if chunk_out.summary is None:
        outs.summary = None
        outs.inferred_length_distribution = None
    else:
        shutil.copy(chunk_out.summary, outs.summary)
        shutil.copy(chunk_out.inferred_length_distribution, outs.inferred_length_distribution)


def main(args, outs):

    stats = main_report_length_mass(args, outs)

    # Write fragment length hist to json
    with open(outs.summary, 'w') as summary_file:
        tenkit.safe_json.dump_numpy(stats, summary_file)

    with open(outs.inferred_length_distribution, 'w') as ld_file:
        tenkit.safe_json.dump_numpy(stats['length_distribution'], ld_file)


def empirical_length_distribution(fragment_df):
    # Generate a fragment length histogram
    len_bin = np.floor_divide(fragment_df.est_len.values, FRAG_LEN_HIST_BIN_SIZE).astype(int) * FRAG_LEN_HIST_BIN_SIZE
    len_bin = len_bin[(fragment_df.num_reads.values > 1).astype('bool')]

    len_bin_counts = {}
    for lb in len_bin:
        len_bin_counts[lb] = len_bin_counts.get(lb, 0) + 1

    len_hist = {str(k):v for (k,v) in len_bin_counts.iteritems()}
    return len_hist

def main_report_length_mass(args, outs):
    tmp_dir = os.path.dirname(outs.summary)

    empty_stats = {
        'alpha': [],
        'alpha_mean': None,
        'alpha_cv': None,
        'mean_frags': None,
        'total_frags': [],
        'length_distribution': {},
        'empirical_length_distribution': {},
        'inferred_mean_length': None,
        'inferred_lw_mean_length': None,
        'inferred_total_mass_ng': None,
        'inferred_bp_per_bc': [],
        'mean_bp_per_bc': 0,
        'occupied_bcs': 0,
        'inferred_number_gems': 0,
    }


    if args.barcodes is None or args.barcode_whitelist is None or not os.path.exists(args.barcodes):
        return empty_stats

    barcode_whitelist = tk_seq.load_barcode_whitelist(args.barcode_whitelist)

    if len(barcode_whitelist) < 1000:
        return empty_stats

    if args.targets_file is None:
        targeted = False
        num_frags = NUM_FRAGS
    else:
        targeted = True
        num_frags = NUM_FRAGS_TARGETED

    bc_df = tenkit.hdf5.read_data_frame(args.barcodes)
    frag_df = tenkit.hdf5.read_data_frame(args.fragments, ['bc', 'chrom', 'start_pos', 'obs_len', 'num_reads', 'est_len'])
    input_num_frags = len(frag_df)

    gem_group = [int(bc.split('-')[1]) for bc in bc_df.bc]
    num_gem_groups = len(set(gem_group))

    # Start with data about all barcodes.
    # First filter out any barcodes that don't have at least 1 molecule that has > 1 read
    # This eliminates most of the background contamination of barcodes
    bc_df = bc_df[bc_df.bc_mean_reads_per_fragment > 1.0].copy()
    bc_df.sort('bc_num_reads', inplace=True)

    # Subset set to the N99 barcodes. (i.e. barcode that account for 99% of reads), and have at least 1 valid fragment
    # A valid fragment must have >= 1 MAPQ30 read and at least 1
    bc_df['cum_reads'] = np.cumsum(bc_df.bc_num_reads)
    prod_bc_thresh = 0.01 * bc_df.bc_num_reads.sum()
    occupied_bcs_df = bc_df[np.logical_and(bc_df.cum_reads > prod_bc_thresh, bc_df.bc_num_fragments > 0)]

    if len(occupied_bcs_df) == 0:
        martian.log_info("No valid barcodes for length/mass inference -- exiting")
        return empty_stats

    # Figure out the subset of BCs likely to be singleton BCs
    # Only run estimation on that subset
    # Infer the expected total GEM count that should have been present
    occupied_bcs = len(occupied_bcs_df)
    total_diversity = len(barcode_whitelist) * num_gem_groups


    # Poisson correction -- we know how many barcodes have >= 1 GEM, and we know
    # how many total barcodes are possible. he Poission distribution to back-calculate
    # The number of GEMs that must have been present.
    # For Chromium there are 4.2M barcodes.
    p_occupied = float(occupied_bcs) / total_diversity
    mean_gems_per_bc = min(100, -np.log(1 - p_occupied))
    p_singleton = scipy.stats.poisson.pmf(1, mean_gems_per_bc)
    n_singleton = p_singleton * total_diversity

    # n_gems gets reported out as 'Gems Detected' in Loupe
    n_gems = int(round(mean_gems_per_bc * total_diversity))


    # Only use the bottom 90% of singleton BCs, to avoid contamination at high end
    bc_df_frags = occupied_bcs_df.sort('bc_num_fragments')
    singleton_bcs = bc_df_frags[int(round(n_singleton*0.1)):int(round(n_singleton*0.9))]

    martian.log_info("Read Count Threshold for Occupied Barcodes: %f" % occupied_bcs_df.iloc[0].bc_num_reads)
    martian.log_info("Occupied Barcodes: %d" % occupied_bcs)
    martian.log_info("Singleton Barcodes: %f" % n_singleton)
    martian.log_info("Number of GEMs in slice used for inference: %d" % len(singleton_bcs))
    martian.log_info("Inferred Number of GEMS: %f" % n_gems)

    # Get empirical fragment length distribution
    obs_len = frag_df.obs_len.values

    # It's possible for multi-read fragments to have a size of zero, which
    # causes a vanishing density - set a lower limit
    obs_len = np.maximum(obs_len, 200)

    empirical_dist = empirical_length_distribution(frag_df)

    # Cap the obs_len at a reasonable value, then set the length bins accordingly
    if targeted:
        max_len_adj_factor = 1.6
    else:
        max_len_adj_factor = 1.3

    # select the max length for the fragment length distribution
    max_len = np.int32(np.percentile(obs_len, 99.97) * max_len_adj_factor)
    max_len = np.maximum(max_len, 100000)
    obs_len = np.minimum(obs_len, max_len, dtype=np.int32)
    max_bin = max_len * 1.01
    bin_data = gen_bin_length(NUM_LENGTH_BINS, min_len=500, max_len=max_bin)

    martian.log_info("Fragments trimmed to max length of %d" % max_len)

    # Select a random subset of BCS to work with
    # Fix random seed so that we get repeatable results
    num_bcs = max(MIN_BCS, float(num_frags) / singleton_bcs.bc_num_fragments.mean())
    np.random.seed(0)
    if len(singleton_bcs) > 0:
        sel_bcs = singleton_bcs.irow(np.random.randint(0,len(singleton_bcs), num_bcs)).copy()
    sel_bcs['bc_id'] = np.arange(1,len(sel_bcs)+1)
    sel_frags = frag_df[frag_df.bc.isin(sel_bcs.bc)].copy()
    sel_frags['bc_string'] = sel_frags.bc.astype('string')
    sel_frags.sort(['bc_string'], inplace=True)
    martian.log_info("Usings %d fragments" % len(sel_frags))

    bc_id_lookup = {}
    for (bc, bc_id) in zip(sel_bcs.bc, sel_bcs.bc_id):
        bc_id_lookup[bc] = bc_id
    # Write out the fragment data for stan to consume
    nbcs = len(sel_bcs)

    obs_len = sel_frags.obs_len.values
    # It's possible for multi-read fragments to have a size of zero, which
    # causes a vanishing density - set a lower limit
    obs_len = np.maximum(obs_len, 200)
    # obs_len for single-read fragments is 1000 in the
    # fragment file -- remap to 0
    obs_len[sel_frags.num_reads.values == 1] = 0.0
    obs_len = np.minimum(obs_len, max_len, dtype=np.int32)

    # Data to be passed to stan
    data = {
        # data sizes
        'N': len(sel_frags),
        'BC': nbcs,

        # Per BC stats
        'bc_observed_frags': sel_bcs.bc_num_fragments,

        # Fragment data: bc_id maps fragments to bc, num_reads, and obs_length fragment stats
        'bc_id': [bc_id_lookup[bc] for bc in sel_frags.bc],
        'num_reads': sel_frags.num_reads,
        'obs_length': obs_len,
    }

    # The number of sizes of the length bins
    data.update(bin_data)

    # Add extra data for targeting if neccesary
    if args.targets_file is not None:
        targets = tk_io.get_target_regions_dict(open(args.targets_file))
        fasta = tenkit.reference.open_reference(args.reference_path)
        ctg_sizes = [(name, len(seq)) for (name, seq) in fasta.items()]
        genome_size = float(sum(l for (name,l) in ctg_sizes))

        gb_size = 1024
        ctg_round_sizes = np.array([math.ceil(float(sz)/gb_size) * gb_size for (name, sz) in ctg_sizes])
        ctg_starts = np.cumsum(np.concatenate([[0], ctg_round_sizes[:-1]]))
        ctg_start_series = p.Series(np.array(ctg_starts,dtype=np.int64), index=[name for (name, l) in ctg_sizes])

        targ_cs_ctgs = []
        on_target_bps = {}
        rsum = 0
        for ((ctg, sz), round_sz) in zip(ctg_sizes, ctg_round_sizes):
            targs = np.zeros(round_sz, dtype=np.int32)
            # Mark bases as targeted
            for (s,e) in targets.get(ctg, []):
                targs[s:e] = 1

            for frag_len in data['bin_length']:
                on_target_chrom = np.zeros(round_sz, dtype=np.int8)

                for (s,e) in targets.get(ctg, []):
                    ss = max(0, s - int(frag_len))
                    ee = min(round_sz, e)
                    on_target_chrom[ss:ee] = 1

                # Determine the probability that a fragment w/ a given length will touch an exon
                on_target_bps[frag_len] = on_target_bps.get(frag_len, 0) + on_target_chrom.sum()
                del on_target_chrom

            # Running sum over chromosomes
            targs_cs = np.cumsum(targs) + rsum
            rsum += np.sum(targs)
            targ_cs_bins = targs_cs[::gb_size].copy()
            del targs
            del targs_cs
            targ_cs_ctgs.append(targ_cs_bins)

        total_target_size = sum((e-s) for regs in targets.values() for (s,e) in regs)
        print "Total target size: %d" % total_target_size
        on_target_fracs = { k:float(v)/genome_size for (k,v) in on_target_bps.items() }
        print on_target_fracs

        # STAN will use this to interpolate the target sizes
        cum_target_bins = np.concatenate(targ_cs_ctgs)

        assert(cum_target_bins.shape[0] == int(np.sum(ctg_round_sizes/gb_size)))

        # Get the position of each fragment on the laid-out genome, with the position decimated by 8
        ctg_starts = ctg_start_series[sel_frags.chrom].values
        stan_pos = ((ctg_starts + sel_frags.start_pos) / 8).astype(np.int32)
        sel_frags['stan_pos'] = stan_pos

        print sel_frags.head(20)

        data['pos'] = sel_frags.stan_pos
        data['genome_size'] = genome_size
        data['gb_size'] = gb_size
        data['GB'] = len(cum_target_bins)
        data['cum_target_bases'] = cum_target_bins

    # Write out the stan input data
    input_fn = os.path.join(tmp_dir, "input.R")
    write_stan_input(input_fn, data)

    # Generate initial values for optimization
    ramp = np.linspace(1,0.1, NUM_LENGTH_BINS)
    ramp = ramp / ramp.sum()

    # assume that fragments with 1 read were 2kb when setting initial alpha
    seen_dna = sel_frags.obs_len.sum() + 2000.0 * (sel_frags.num_reads == 1).sum()
    mean_alpha = float(sel_frags.num_reads.sum()) / seen_dna

    frags_mu = sel_bcs.bc_num_fragments.mean()

    # Initial values of parameters to be estimated by Stan
    init_data = {
        # BC amp rate
        'alpha': [mean_alpha] * nbcs,

        # Length distribution
        'theta': list(ramp),

        # Average number of fragments
        'mean_frags': frags_mu,

        # Number of unobserved fragments
        'bc_unobserved_frags': [100] * nbcs,
        'read_disp': 10,
        'amp_length_k': 1.0/200000,
    }

    init_fn = os.path.join(tmp_dir, "init.R")
    write_stan_input(init_fn, init_data)

    # check if we have valid data for stan
    # need some observed fragments, and a minimal reads / fragments
    mean_rpf = sel_frags.num_reads.mean()
    martian.log_info("Mean LPM of molecules selected for inference: %f" % mean_rpf)

    success = 0
    if len(sel_frags) > 0 and mean_rpf > MIN_RPF and (not targeted or total_target_size >= MIN_TARGET_SIZE):
        success = run_model(tmp_dir, targeted)
    else:
        if targeted and total_target_size < MIN_TARGET_SIZE:
            martian.log_info("Target size is too small for length/mass inference: %d" % total_target_size)

        if len(sel_frags) == 0:
            martian.log_info("Aborting length-mass inference: no fragments")

        if mean_rpf < MIN_RPF:
            martian.log_info("Reads per fragment too low for length-mass inference: %f" % mean_rpf)

    if success:
        res = load_stan_output(os.path.join(tmp_dir, "output.csv"))

        # If targeted, adjust the fragment length distribution and mass according to the fragment
        # visibility function
        if targeted:
            theta = res['theta']
            bl = data['bin_length']
            vis_func = np.array([on_target_fracs[l] for l in bl])
            print vis_func
            adj_theta = theta / vis_func
            adj_theta = adj_theta / adj_theta.sum()

            missing_factor = 1.0 / (adj_theta * vis_func).sum()

            # Put back in the adjusted values
            res['theta'] = adj_theta
            res['mean_frags'] = missing_factor * res['mean_frags']
            res['bc_total_frags'] = missing_factor * res['bc_total_frags']


        # print the mass distribution, alpha distributions
        mean_length = (data['bin_length'] * res['theta']).sum()
        mean_length_weighted = np.average(data['bin_length'], weights=data['bin_length']*res['theta'])

        # Mass conversion
        ng_per_bp = 1.025e-12

        bases_per_bc = res['bc_total_frags'] * mean_length
        total_bases = res['bc_total_frags'].mean() * mean_length * n_gems
        total_mass_ng = total_bases * ng_per_bp

        # calculation
        bp_per_ng = 9.76e11

        # try to calc input mass
        #z2_vol_per_gem - ufluidcs number, corrected for empty GEMS
        #bp_per_gem = loaded_mass * bp_per_ng * z2_vol_per_gem / total_z2_vol_input
        # z2_vol_per_gem = 144 pL
        # total_z2_vol_input = 65uL
        # Fixme -- product configuration needs to be passed in & fixed for future products
        fluidics_params = FLUIDICS_PARAMS['Chromium']
        loaded_mass = np.mean(bases_per_bc) * fluidics_params['total_z2_vol_input'] / bp_per_ng / fluidics_params['z2_vol_per_gem']

        # Me: magic number, David: empirically derived correction factor
        DENATURATION_FACTOR = 1.6

        # Ad-hoc correction for the apparent 'denaturation' of the input material, which leads to double counting on input DNA
        corrected_loaded_mass = loaded_mass / DENATURATION_FACTOR

        stats = {
            'alpha': list(res['alpha']),
            'alpha_mean': np.mean(res['alpha']),
            'alpha_cv': tk_stats.robust_divide(np.std(res['alpha']), np.mean(res['alpha'])),
            'mean_frags': res['mean_frags'],
            'total_frags': res['bc_total_frags'],
            'length_distribution': { str(l):frac for (l,frac) in zip(data['bin_length'], input_num_frags*res['theta']) },
            'empirical_length_distribution': empirical_dist,
            'inferred_mean_length': mean_length,
            'inferred_lw_mean_length': mean_length_weighted,
            'inferred_total_mass_ng': total_mass_ng,
            'inferred_bp_per_bc': bases_per_bc,
            'mean_bp_per_bc': np.mean(bases_per_bc),
            'loaded_mass_ng': loaded_mass,
            'corrected_loaded_mass_ng': corrected_loaded_mass,
            }
    else:

        len_dist_default = { str(k):1.0/k for k in data['bin_length'] }

        stats = {
            'alpha': [],
            'alpha_mean': None,
            'alpha_cv': None,
            'mean_frags': None,
            'total_frags': [],
            'length_distribution': len_dist_default,
            'empirical_length_distribution': empirical_dist,
            'inferred_mean_length': None,
            'inferred_lw_mean_length': None,
            'inferred_total_mass_ng': None,
            'inferred_bp_per_bc': [],
            'mean_bp_per_bc': None,
            'loaded_mass_ng': None,
            'corrected_loaded_mass_ng': None,
            }

    stats['occupied_bcs'] = occupied_bcs
    stats['inferred_number_gems'] = n_gems
    return stats



def run_model(dirname, targeted = False):
    if targeted:
        exe = "len_mass_model_targeted"
        iters = 600
    else:
        exe = 'len_mass_model'
        iters = 900

    args = [ exe, 'optimize',
            'iter=%d' % iters,
            'data',
            "file=" + os.path.join(dirname, "input.R"),

            'output',
            'file=' + os.path.join(dirname, "output.csv"),
            "refresh=10",
            "init=" + os.path.join(dirname, "init.R"),
            "random",  "seed=4072180573",
            ]

    print " ".join(args)
    proc = subprocess.Popen(args, stdout=subprocess.PIPE)
    (stdoutdata, stderrdata) = proc.communicate()
    martian.log_info(stdoutdata)
    martian.log_info(stderrdata)
    return (proc.returncode == 0)
