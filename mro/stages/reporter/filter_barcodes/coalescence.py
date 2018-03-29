#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import numpy as np
import scipy
import scipy.stats
import scipy.sparse
import tenkit.hdf5 as kt_hdf
import tenkit.bam as tk_bam
import sys
import tenkit.hdf5
import cProfile
import tenkit.pandas as p
import martian
import union_find

GENOME_BINS=60000 # = roughly 50kb bin length for human

class BcSimilarity:

    def __init__(self, fragments_fn, bcs_to_use, bam_fn):
        ''' Compute sparse barcode x genome bin coverage matrix.
            Each row is normalized to 1, so that the expected
            overlap for uncorrelated barcodes is 1 '''

        self.uf = union_find.UnionFind()

        bam_in = tk_bam.create_bam_infile(bam_fn)

        # Choose the set of fragments to use
        if type(fragments_fn) is str or type(fragments_fn) is unicode:
            fragments = self.load_fragments_filtered(fragments_fn, bcs_to_use)
        elif type(fragments_fn) is p.DataFrame:
            fragments = fragments_fn
        else:
            raise Exception("unrecognized fragments_fn argument type: %s, must be filename or pandas.DataFrame" % str(type(fragments_fn)))

        # Setup genome bins
        genome_length = sum(bam_in.lengths)
        bin_size = max(1, genome_length/GENOME_BINS)
        chrom_bins = np.ceil(np.array([float(l) / bin_size for l in bam_in.lengths]))
        total_bins = chrom_bins.sum()
        start_bin = np.concatenate([[0], np.cumsum(chrom_bins)[:-1]])
        chrom_map = {c:idx for (idx,c) in enumerate(bam_in.references)}

        npartitions = len(bcs_to_use)

        # Number the selected barcodes -- the assigned number is their row in the BC-bin matrix
        bcs = fragments.bc.values
        bc_ids = {}
        self.bcs_to_use = []
        c = 0
        for bc in bcs:
            if bc_ids.has_key(bc):
                continue
            else:
                self.bcs_to_use.append(bc)
                bc_ids[bc] = c
                c += 1

        martian.log_info("making sparse matrix")

        indexes = np.empty((2,len(fragments)), dtype=np.int32)
        data = np.ones((len(fragments),), dtype=np.float32)

        chroms = fragments.chrom.values
        pos_bin = fragments.start_pos.values / bin_size

        for fidx in range(len(fragments)):
            chrom_id = chrom_map[chroms[fidx]]
            which_bin = start_bin[chrom_id] + pos_bin[fidx]
            which_bc = bc_ids[bcs[fidx]]
            indexes[0,fidx] = which_bc
            indexes[1,fidx] = which_bin

        mat = scipy.sparse.csr_matrix((data, indexes), shape=(npartitions, total_bins), dtype=np.float32)
        # If there are multiple fragments for the same BC in the same bin, the csr_matrix constructor above will sum them up, leading
        # to entries greater than 1.  Cap everything at 1.
        mat.data = np.ones(mat.data.shape, dtype=mat.dtype)

        '''
        mat1 = scipy.sparse.lil_matrix((npartitions, total_bins), dtype=np.float32)
        bc_grps = fragments.groupby(["bc"])
        bc_count = 0

        # For each barcode, mark the genome bins covered by a fragment
        for (bc, bc_grp) in bc_grps:
            # Track the reads per fragment in tested partitions for reporting
            l = len(bc_grp)
            bins = np.zeros(l, dtype=np.int32)
            chroms = bc_grp.chrom.values
            starts = bc_grp.start_pos.values
            pos_bin = starts / bin_size

            for i in range(l):
                chrom_id = chrom_map[chroms[i]]
                which_bin = start_bin[chrom_id] + pos_bin[i]
                bins[i] = which_bin

            mat1[bc_count, bins] = 1.0
            bc_count += 1

            if bc_count % 1000 == 0:
                print bc_count
        '''
        eps = 0.0001

        # Get the genome bin occupancy
        genome_bin_counts = np.array((mat > np.float32(0)).sum(axis=0)).flatten().astype('float') # total BC counts per bin
        high_cov_threshold = np.percentile(genome_bin_counts, 99.5)

        # switch off high-coverage bins -- set them to eps (a small nonzero number so we can distinguish them)
        high_cov_bins = np.where(genome_bin_counts > high_cov_threshold)[0]
        (r,c) = mat.nonzero()
        martian.log_info("removing %d bins" % len(high_cov_bins))
        for hc_bin in high_cov_bins:
            mat.data[c == hc_bin] = eps

        # Recalculate the genome bins distribution
        genome_bin_counts = np.array((mat > (2*eps)).sum(axis=0)).flatten().astype('float')
        martian.log_info("Genome Bin Coverage  mean: %f  99.95th percentile: %f" % (genome_bin_counts.mean(), high_cov_threshold))

        # Adjust for 'effective genome size' based on the distribution over bins
        # i.e. more skewed distribution -> fewer effective bins
        effective_bins_factor = ((genome_bin_counts / genome_bin_counts.sum())**2).sum()
        self.effective_genome_bins = 1.0/effective_bins_factor
        martian.log_info("Effective Number of Genome Bins = %f" % self.effective_genome_bins)

        self.mat = mat
        martian.log_info("done __init__")


    def load_fragments_filtered(self, fn, bcs_to_use):
        ''' Load fragment data for coalescence calculation '''

        martian.log_info("loading fragment data")

        def fragment_filter(frags):
            #return np.logical_and(frags.num_reads > 1, frags.bc.isin(bcs_to_use))
            return frags.bc.isin(bcs_to_use)

        frags = kt_hdf.read_data_frame_filtered(fn, fragment_filter, query_cols=['bc', 'num_reads', 'chrom', 'start_pos'])
        return frags


    def coalescence_analysis(self, min_cluster_size = 2, fpr = 1.0):
        ''' Compute the BC-BC overlap matrix, threshold it and convert to a graph, and report large clusters '''
        # Sizes
        (num_bcs, total_bins) = self.mat.shape

        # number of fragment on each bc
        self.num_frags = np.array((self.mat > 0.0).sum(axis=1)).flatten()
        martian.log_info("mean num frags: %f,  median: %f" % (self.num_frags.mean(), np.median(self.num_frags)))

        # compute an initial threshold for coalescence
        # first, compute the expected number of overlaps between two BCs having the mean number of fragments
        expected_overlaps = (self.num_frags.mean()**2) / self.effective_genome_bins
        # now use a ~5 sigma threshold to represent an initial cutoff for significance
        # note: this is more informative for GemCode than Chromium, because there are many more fragments per BC
        # and thus the expected number of overlaps due to chance is much higher. for Chromium, the threshold will usually be 2.
        overlap_threshold = np.float32(max(2, round(expected_overlaps + 5 * np.sqrt(expected_overlaps))))
        martian.log_info("expected overlaps: %f -- using overlap threshold: %f" % (expected_overlaps, overlap_threshold))

        # Chunk out matrix in x and y and find significant overlaps in each pair of chunks
        bc_bin_size = 1000
        bc_bins = np.arange(0, num_bcs, bc_bin_size)

        # Choose a p-value that accounts for many comparisons
        n_comparisons = num_bcs**2
        pvalue_cut = fpr / n_comparisons
        martian.log_info("using pvalue cutoff: %e" % pvalue_cut)

        for x in bc_bins:
 
            martian.log_info("BC number: %d on x axis" % x)
            for y in bc_bins:
                # calculation is symmetric -- don't do below the diagonal
                if y < x:
                    continue

                self.window_intersection_slices(x, y, bc_bin_size, pvalue_cut, overlap_threshold)


        # Delete the fragment matrix to save memory
        del self.mat
        self.mat = None

        martian.log_info("Finding connected components")
        clusters = self.clusters()
        bad_bcs = []
        cluster_idxs = []
        cluster_size = []

        martian.log_info("Making bad BC DataFrame")
        for (cluster_idx, (_, nodes)) in enumerate(clusters.iteritems()):
            if len(nodes) < min_cluster_size:
                continue
            bad_bcs.extend(nodes)
            cluster_idxs.extend([cluster_idx] * len(nodes))
            cluster_size.extend([len(nodes)] * len(nodes))

        bad_bc_tbl = p.DataFrame({'bc': bad_bcs, 'cluster_id': cluster_idxs, 'cluster_size': cluster_size})
        return bad_bc_tbl

    def window_intersection_slices(self, bc_offset1, bc_offset2, bc_bin_size, pvalue_threshold, overlap_threshold):
        '''Count the window intersection for bc slices of the full matrix'''

        num_bcs, num_windows = self.mat.shape
        slice1 = self.mat[bc_offset1:min(num_bcs, bc_offset1 + bc_bin_size)]
        slice2 = self.mat[bc_offset2:min(num_bcs, bc_offset2 + bc_bin_size)].transpose()

        # convert lil matrix to csr and multiply by it's transpose
        shared_fragments = slice1 * slice2

        # overlap matrix is normalized so that the average value for uncorrelated BCs
        # is 1.  Use that to pick a threshold of 'interesting' BC pairs

        # Take a very conservative bite to reduce false positives
        overlapping_bcs  = shared_fragments > overlap_threshold
        _i1,_i2 = overlapping_bcs.nonzero()
        first = _i1 < _i2

        # scipy.sparse has bug for empty selection -- bail out here to avoid it
        if first.sum() == 0:
            return []

        shared_fragments_pass = np.array(shared_fragments[_i1[first], _i2[first]]).flatten()

        i1s = _i1[first] + bc_offset1
        i2s = _i2[first] + bc_offset2

        for (elem_idx, (i1, i2)) in enumerate(zip(i1s, i2s)):

            if self.same_cluster(i1, i2):
                continue

            # Compute the p-value of the overlaps under a null model with uncorrelated fragments
            lmbda = float(self.num_frags[i1] * self.num_frags[i2]) / self.effective_genome_bins
            overlaps = shared_fragments_pass[elem_idx]

            # Do a quick check that we have a significant hit before doing detailed p-value calculation
            if overlaps < lmbda + 6 * np.sqrt(lmbda):
                continue

            pvalue = scipy.stats.poisson.sf(overlaps, lmbda)
            if pvalue < pvalue_threshold:
                self.detected_overlap(self.bcs_to_use[i1], self.bcs_to_use[i2])

    def same_cluster(self, i1, i2):
        if self.uf.has_item(i1) and self.uf.has_item(i2):
            return self.uf[i1] == self.uf[i1]

    def detected_overlap(self, i1, i2):
        self.uf.union(i1, i2)

    def clusters(self):
        clusters = {}
        for item in self.uf:
            cluster = self.uf[item]
            node_list = clusters.setdefault(cluster, [])
            node_list.append(item)

        return clusters

def go():
    path = sys.argv[1]
    fragment_fn = path + "/fork0/files/fragments.h5"
    barcode_fn = path + "/fork0/files/barcodes.h5"

    bc = tenkit.hdf5.read_data_frame(barcode_fn)

    gbcs = bc[np.logical_and(bc.bc_mean_reads_per_fragment > 1.2, bc.bc_num_reads > 50)]
    print "number of good bcs: %d" % len(gbcs)
    bcs_to_use = list(gbcs.bc)
    bam_fn = sys.argv[2]

    bcSim = BcSimilarity(fragment_fn, bcs_to_use, bam_fn)
    (bad_bcs_tbl, graph) = bcSim.coalescence_analysis()

    import cPickle
    with open("coa.pck", "wb") as f:
        cPickle.dump(graph, f, protocol=2)

    bad_bcs_tbl.to_csv("coa.tsv", sep="\t")

# Script entry point for debugging / exploratory analysis
if __name__ == '__main__':

    profile = cProfile.Profile()
    profile.enable()
    try:
        profile.run("go()")
    finally:
        profile.disable()
        profile.dump_stats("profile.out")


