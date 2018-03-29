#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Attach phasing information to BAM records
#
import subprocess
import tenkit.tabix as tk_tabix
import tenkit.hdf5 as tk_hdf5
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
from tenkit.constants import FRAGMENT_PHASING_THRESHOLD, PHASE_SET_BAM_TAG, HAPLOTYPE_BAM_TAG, MOLECULE_ID_BAM_TAG, PHASING_CONF_BAM_TAG
import martian
import tenkit.stats as tk_stats
import json
import tenkit.safe_json

__MRO__ = """
stage ATTACH_PHASING(
    in  bam     input,
    in  tsv.gz  fragment_phasing,
    in  h5      fragments,
    out bam     phased_possorted_bam,
    out bam.bai phased_possorted_bam_index,
    out int     total_reads,
    out int     phased_reads,
    out int     molecule_tagged_reads,
    src py      "stages/reads/attach_phasing",
) split using (
    in  string  chunk_start,
    in  string  chunk_end,
)
"""

def split(args):
    bam_in = tk_bam.create_bam_infile(args.input)
    chunk_defs = tk_bam.chunk_bam_records(bam_in, chunk_bound_key=None)
    for c in chunk_defs:
        c["__mem_gb"] = 4
    return {'chunks': chunk_defs}

WINDOW_SIZE=500000

def main(args, outs):
    chunk_start = args.chunk_start
    chunk_end = args.chunk_end

    bam_in = tk_bam.create_bam_infile(args.input)
    reads = tk_bam.read_bam_chunk(bam_in, (chunk_start, chunk_end))
    pgs = [
        tk_bam.make_pg_header(martian.get_pipelines_version(), "attach_phasing"),
        tk_bam.make_terminal_pg_header(martian.get_pipelines_version())
    ]
    # dont duplicate header if already there, this is for developer testing purposes
    PG = bam_in.header['PG']
    for item in PG:
        if item['ID'] == "attach_phasing":
            pgs = []
    bam_out, _ = tk_bam.create_bam_outfile(outs.phased_possorted_bam, None, None, template=bam_in, pgs=pgs)

    # File with contig phasing information
    if args.fragment_phasing is not None:
        frag_phasing = tk_tabix.create_tabix_infile(args.fragment_phasing)
    else:
        frag_phasing = None

    if args.fragments is not None:
        # Fragments file for global molecule id
        frag_id_reader = tk_hdf5.DataFrameReader(args.fragments)
    else:
        frag_id_reader = None

    # Phasing data
    ph_db = None
    ph_db_chrom = None
    ph_db_start = None
    ph_db_end = None

    # Fragment data - for global molecule id
    fr_db = None
    fr_db_chrom = None
    fr_db_start = None
    fr_db_end = None

    total_reads = 0
    phased_reads = 0
    molecule_tagged_reads = 0

    for r in reads:
        chrom = bam_in.references[r.tid]
        pos = r.pos
        bc = tk_io.get_read_barcode(r)

        total_reads += 1
        tags = r.tags

        # Strip out RX and QX tags
        #strip_tags = [RAW_BARCODE_TAG, RAW_BARCODE_QUAL_TAG]
        # Actually don't strip
        strip_tags = []
        tags = [(tag, value) for (tag, value) in tags if (tag not in strip_tags)]

        # Fetch from the fragments file to get records that should cover many future reads
        # fragment phasing file may not exist in ALIGNER only pipeline - may need to skip
        if frag_phasing is not None:
            if ph_db is None or chrom != ph_db_chrom or pos < ph_db_start or pos > ph_db_end:
                ph_db, (ph_db_chrom, ph_db_start, ph_db_end) = get_frag_phasing_db(frag_phasing, chrom, pos, window=WINDOW_SIZE)

            if bc is not None and ph_db.has_key(bc):
                frags = ph_db[bc]
                # See if we having phasing for this fragment
                valid_phasing = [x for x in frags if x['start'] <= r.pos and x['end'] > r.pos]
                assert(len(valid_phasing) < 2)
                if len(valid_phasing) == 1:
                    phased_reads += 1
                    read_phasing = valid_phasing[0]
                    tags.append((PHASE_SET_BAM_TAG, read_phasing['ps']))
                    tags.append((HAPLOTYPE_BAM_TAG, read_phasing['hap']))
                    tags.append((PHASING_CONF_BAM_TAG, read_phasing['pc']))

        if frag_id_reader is not None:
            # Fetch from the fragments file to get records that should cover many future reads
            if fr_db is None or chrom != fr_db_chrom or pos < fr_db_start or pos > fr_db_end:
                fr_db, (fr_db_chrom, fr_db_start, fr_db_end) = get_molecule_id_db(frag_id_reader, chrom, pos, window=WINDOW_SIZE)

            if bc is not None and fr_db.has_key(bc):
                frags = fr_db[bc]
                # See if we having phasing for this fragment
                molecule_ids = [x for x in frags if x['start'] <= r.pos and x['end'] > r.pos]
                assert(len(molecule_ids) < 2)
                if len(molecule_ids) == 1:
                    molecule_tagged_reads += 1
                    molecule_id = molecule_ids[0]
                    tags.append((MOLECULE_ID_BAM_TAG, molecule_id['molecule_id']))


        r.tags = tags
        bam_out.write(r)

    bam_out.close()
    outs.total_reads = total_reads
    outs.phased_reads = phased_reads
    outs.molecule_tagged_reads = molecule_tagged_reads

def get_frag_phasing_db(fragment_file, chrom, pos, window=int(1e6)):
    db = {}

    try:
        frag_phase_iter = fragment_file.fetch(chrom, pos, pos + window)
    except IndexError:
        frag_phase_iter = []
    except KeyError:
        frag_phase_iter = []
    except ValueError:
        frag_phase_iter = []


    for frag_line in frag_phase_iter:
        frag = frag_line.strip().split('\t')
        start, end = int(frag[1]), int(frag[2]) # contig coordinates
        ps, ps_start, ps_end = int(frag[3]), int(frag[4]), int(frag[5])
        bc = frag[6]
        p_hap1 = float(frag[7])
        p_hap2 = float(frag[8])

        if p_hap1 > FRAGMENT_PHASING_THRESHOLD or p_hap2 > FRAGMENT_PHASING_THRESHOLD:
            if p_hap1 > p_hap2:
                hap = 1
                pc = tk_stats.qual_from_prob_correct(p_hap1)
            else:
                hap = 2
                pc = tk_stats.qual_from_prob_correct(p_hap2)

            phase_set_frag = {
                'start': max(start, ps_start),
                'end': min(end, ps_end),
                'ps': ps,
                'hap': hap,
                'pc': pc
            }

            bc_list = db.setdefault(bc, [])
            bc_list.append(phase_set_frag)

    return db, (chrom, pos, pos + window)

def get_molecule_id_db(fragment_reader, chrom, pos, window=int(1e6)):
    db = {}

    frags = fragment_reader.query((chrom, pos, pos + window), query_cols=['bc', 'start_pos', 'end_pos', 'molecule_id'])

    for k,r in frags.iterrows():

        molecule_id = {
                'start': r.start_pos,
                'end': r.end_pos,
                'molecule_id': r.molecule_id,
            }

        bc_list = db.setdefault(r.bc, [])
        bc_list.append(molecule_id)

    return db, (chrom, pos, pos + window)



def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    # Concatenate chunks
    if len(chunk_outs) == 1:
        subprocess.call(['mv', chunk_outs[0].phased_possorted_bam, outs.phased_possorted_bam])
    else:
        tk_bam.concatenate(outs.phased_possorted_bam, [out.phased_possorted_bam for out in chunk_outs])
    tk_bam.index(outs.phased_possorted_bam)
    outs.phased_possorted_bam_index = outs.phased_possorted_bam + ".bai"

    total_reads = 0
    phased_reads = 0
    molecule_tagged_reads = 0
    for chunk_out in chunk_outs:
        total_reads += chunk_out.total_reads
        phased_reads += chunk_out.phased_reads
        molecule_tagged_reads += chunk_out.molecule_tagged_reads

    outs.total_reads = total_reads
    outs.phased_reads = phased_reads
    outs.molecule_tagged_reads = molecule_tagged_reads

    fract_reads_phased = tk_stats.robust_divide(float(phased_reads), float(total_reads))
    fract_reads_molecule_id = tk_stats.robust_divide(float(molecule_tagged_reads), float(total_reads))

    stats = {
        "fract_reads_phased": fract_reads_phased,
        "fract_reads_molecule_id": fract_reads_molecule_id,
        }

    with open(outs.summary, 'w') as summary_file:
        json.dump(tenkit.safe_json.json_sanitize(stats), summary_file)
