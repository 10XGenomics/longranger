#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import re
import socket
import subprocess
import vcf
import martian
import tenkit.preflight as tk_preflight
import tenkit.bio_io as tk_io
import tenkit.reference
import traceback

__MRO__ = """
stage PHASER_SVCALLER_PREFLIGHT(
    in string sample_id,
    in map[]  sample_def,
    in string sex,
    in bed    targets,
    in bed    target_blacklist,
    in string restrict_locus,
    in string vc_mode,
    in string reference_path,
    in vcf    vc_ground_truth,
    in int    sv_min_qv,
    in bedpe  sv_ground_truth,
    in bool   check_executables,
    src py     "stages/preflight/phaser_svcaller",
)
"""

def main(args, outs):
    hostname = socket.gethostname()

    # Sample ID / pipestance name
    if args.sample_id is not None:
        if not re.match("^[\w-]+$", args.sample_id):
            martian.exit("Sample name may only contain letters, numbers, underscores, and dashes: " + args.sample_id)

    # FASTQ input
    for sample_def in args.sample_def:
        #if not tk_preflight.check_is_chromium(sample_def):
        #    martian.exit("This version of Longranger does not support GemCode data. Please use Longranger 1.2 instead.")
        read_path = sample_def["read_path"]
        if not read_path:
            martian.exit("Must specify a read_path containing FASTQs.")
        if not read_path.startswith('/'):
            martian.exit("Specified FASTQ folder must be an absolute path: %s" % read_path)
        if not os.path.exists(read_path):
            martian.exit("On machine: %s, specified FASTQ folder does not exist: %s" % (hostname, read_path))
        if not os.access(read_path, os.X_OK):
            martian.exit("On machine: %s, longranger does not have permission to open FASTQ folder: %s" % (hostname, read_path))
        if not os.listdir(read_path):
            martian.exit("Specified FASTQ folder is empty: " + read_path)

        library_id = sample_def.get("library_id")
        if library_id is not None:
            if not re.match("^[\w-]+$", library_id):
                martian.exit("Library name may only contain letters, numbers, underscores, and dashes: " + library_id)

        lanes = sample_def["lanes"]
        if lanes is not None:
            for lane in lanes:
                if not is_int(lane):
                    martian.exit("Lanes must be a comma-separated list of numbers.")

        ok, msg = tk_preflight.check_sample_indices(sample_def)
        if not ok:
            martian.exit(msg)

    # Reference
    MAX_CONTIGS = 1000
    ok, msg = tk_preflight.check_refdata(args.reference_path, MAX_CONTIGS)
    if ok:
        martian.log_info(msg)
    else:
        martian.exit(msg)

    # Sex (given reference)
    if args.sex is not None:
        if args.sex.lower() not in ["m", "male", "f", "female"]:
            martian.exit("Sex of sample must be 'm', 'male', 'f', or 'female'.")
    else:
        if tenkit.reference.load_male_chromosomes(args.reference_path) == None:
            martian.exit("Must specify sex of sample, or use a reference package that includes a sex_chromosomes.tsv file.\nFor more details, see http://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/references")

        ref = tenkit.reference.open_reference(args.reference_path)
        male_chrom = tenkit.reference.load_male_chromosomes(args.reference_path)
        for m in male_chrom:
            if m not in ref:
                martian.exit("Reference issue in sex_chromosomes.tsv. Male-specific chromosome '%s' does not exist in reference" % m)

        auto_chrom = tenkit.reference.load_autosomal_chromosomes(args.reference_path)
        if auto_chrom is None:
            martian.exit("No autosomal chromosome listed in sex_chromosomes.tsv. Please list an autosomal chromosome to use as a reference for sex determination")

        for a in auto_chrom:
            if a not in ref:
                martian.exit("Reference issue in sex_chromosomes.tsv. Autosomal chromosome '%s' does not exist in reference" % a) 

    # Open file handles limit - per LONGRANGER-1758, only check this on the execution machine.
    # We can tell if we're on the execution machine by looking at args.check_executables
    if args.check_executables:
        ok, msg = tk_preflight.check_open_fh()
        if not ok:
            martian.exit(msg)

    # Targets
    if args.targets is not None:
        tk_preflight.check_file("targets", args.targets, hostname)
        tk_preflight.check_bed(args.targets, args.reference_path)

        if args.target_blacklist is None:
            print "\nWARNING: You selected targeted mode but did not provide a --cnvfilter.\nPlease note this may result in a high number of false positive CNV calls.\nFor more details, see http://support.10xgenomics.com/genome-exome/software\n"

    # Target blacklist
    if args.target_blacklist is not None:
        tk_preflight.check_file("cnvfilter", args.target_blacklist, hostname)
        tk_preflight.check_bed(args.target_blacklist, args.reference_path)

    # Restrict locus
    if tenkit.reference.is_tenx(args.reference_path):
        if args.restrict_locus is not None:
            if not re.match("^chr[A-Za-z0-9]{1,2}:[0-9]+\.\.[0-9]+$", args.restrict_locus):
                martian.exit("restrict_locus must be of the form 'chrXX:start..end'.")

    # Pre-called
    if args.vc_precalled is not None:
        tk_preflight.check_file("pre-called VCF", args.vc_precalled, hostname)
        check_vcf(args.vc_precalled, args)

    # VC mode
    if not re.match("^(disable|freebayes|gatk:/.*\.jar|precalled:/.*\.vcf)$", args.vc_mode):
        martian.exit("vc_mode must be of the form 'freebayes', 'gatk:/path/to/gatk_jar_file.jar', 'disable'.")

    if args.vc_precalled is None and args.vc_mode == "disable":
        martian.exit("Because you have not provided a pre-called VCF file, variant calling cannot be disabled.")

    vc_args = args.vc_mode.split(":")
    vc_mode = vc_args[0]
    if vc_mode == "precalled":
        if args.vc_precalled is not None:
            martian.exit("Please specify a pre-called VCF file using only one method.")
        precalled_vars_path = vc_args[1]
        tk_preflight.check_file("pre-called VCF", precalled_vars_path, hostname)
        check_vcf(precalled_vars_path, args)
    elif vc_mode == "gatk":
        jar_path = vc_args[1]
        if not jar_path.startswith('/'):
            martian.exit("Specified GATK jar file must be an absolute path: %s" % jar_path)
        if not os.path.exists(jar_path):
            martian.exit("On machine: %s, specified GATK jar file does not exist: %s" % (hostname, jar_path))
        if os.path.isdir(jar_path):
            martian.exit("Please specify a GATK jar file, not a folder.")
        if args.check_executables:
            check_gatk(jar_path, hostname)

        check_gatk_ref(args.reference_path)

    # VC ground truth
    if args.vc_ground_truth is not None:
        tk_preflight.check_file("VCF ground truth", args.vc_ground_truth, hostname)
        check_vcf(args.vc_ground_truth, args)

    # SV min QV
    if args.sv_min_qv is not None and args.sv_min_qv < 0:
        martian.exit("sv_min_qv must be a positive integer.")

    # SV ground truth
    if args.sv_ground_truth is not None:
        tk_preflight.check_file("SV ground truth", args.sv_ground_truth, hostname)

    martian.log_info(tk_preflight.record_package_versions())


def check_vcf(filename, args):
    fasta = tenkit.reference.open_reference(args.reference_path)
    record_cap = 1000
    record_cursor = 0
    lines = 0
    with open(filename, 'r') as vcf_file:
        for line in vcf_file:
            if lines == 0 and (not line.startswith("##fileformat=VCFv4.")):
                martian.exit(filename + " does not have a proper header. First line should begin with ##fileformat=VCFv4.")
            if not line.startswith("#"):
                break
            lines += 1

    with open(filename, 'r') as f:
        try:
            vcf_iter = vcf.Reader(f)
        except:
            trace = traceback.format_exc()
            martian.exit(filename+" failed on parsing with PyVCF. Traceback:\n"+trace)

        while True:
            try:
                record = vcf_iter.next()
            except StopIteration:
                break
            except:
                trace = traceback.format_exc()
                martian.exit(filename+" failed on parsing with PyVCF. Approximate line number of failure occured at "+str(lines+record_cursor+1)+". Traceback:\n"+trace)
            try:
                record_str = "\nErrored on record " + str(record) + " in file " + filename + " Approximate line number "+str(lines+record_cursor+1)
            except:
                martian.exit(filename+" failed on parsing with pyvcf at approximate line number "+str(lines+record_cursor+1)+". Traceback:\n"+traceback.format_exc())

            # Check for multiple sample columns.
            if len(record.samples) != 1:
                martian.exit("The supplied VCF file contains multiple samples, which is not currently supported: " + str(record.samples))
            try:
                chrom = tk_io.get_record_chrom(record)
                ref = tk_io.get_record_ref(record)
                alt_alleles = tk_io.get_record_alt_alleles(record)
            except:
                martian.exit(filename+" failed on parsing with pyvcf at approximate line number "+str(lines+record_cursor+1)+". Traceback:\n"+traceback.format_exc())

            # Check for chromosome name that doesn't start with 'chr'.
            if tenkit.reference.is_tenx(args.reference_path) and tenkit.reference.get_genome(args.reference_path) == "10X_hg19_ucsc":
                if not chrom.startswith('chr'):
                    martian.exit("The supplied VCF file does not use UCSC-style 'chrX' chromosome names, and this is not currently supported."+record_str)

            # Check that chromosome exists in reference
            if not chrom in fasta:
                martian.exit("The supplied VCF file contains chromosomes not found in reference genome."+record_str)

            # Check that ref allele exists
            if ref is not None:
                if ref == "." or ref == "":
                    martian.exit("The supplied VCF file contains entries with . or missing reference alleles."+record_str)
            else:
                martian.exit("The supplied VCF file contains entries with missing reference alleles."+record_str)

            # Check ref allele is upper case
            if ref != ref.upper():
                martian.exit("The supplied VCF file contains entries with lower case or mixed case reference alleles."+record_str)

            # Check that alt allele isnt empty or '.'
            if alt_alleles is not None:
                for allele in alt_alleles:
                    if allele is None or allele == '.':
                        martian.exit("The supplied VCF file contains entries where ALT allele is either empty or '.'"+record_str)
                    elif allele != allele.upper():
                        martian.exit("The supplied VCF file contains entries with lower case or mixed case alleles."+record_str)
            else:
                martian.exit("The supplied VCF file contains entries with no ALT alleles." + record_str)

            record_cursor += 1
            if record_cursor >= record_cap:
                break
    with open("temp.vcf",'w') as temp:
        subprocess.check_call(['head','-n','3000',filename],stdout=temp)
    with open("temp2.vcf",'w') as temp2:
        try:
            subprocess.check_call(['vcfallelicprimitives','--keep-info','-t','VCFALLELICPRIMITIVE','temp.vcf'], stdout = temp2)
        except:
            trace = traceback.format_exc()
            martian.exit(filename+" failed on parsing with vcfallelicprimitives. Traceback:\n"+trace)
    with open(os.devnull, "w") as fnull:
        try:
            subprocess.check_call(['bcftools', 'filter', 'temp2.vcf'],stdout=fnull)
        except:
            trace = traceback.format_exc()
            martian.exit(filename+" failed on parsing with vcfallelicprimitives or bcftools. Traceback:\n"+trace)
    subprocess.check_call(['rm', 'temp.vcf', 'temp2.vcf'])


def is_int(s):
    try:
        int(s)
    except ValueError:
        return False
    return True


def check_gatk(jar_path, hostname):
    try:
        subprocess.check_call(["which", "java"])
    except subprocess.CalledProcessError:
        martian.exit("Java executable not found on PATH on %s." % hostname)
    try:
        output = subprocess.check_output(["java", "-jar", jar_path, "--version"])
        if "3.3" not in output and "3.4" not in output and "3.5" not in output and "3.7" not in output and "3.8" not in output and "3.9" not in output:
            martian.exit("Version of GATK must be 3.3-3.8, or 4 except 3.6. Version detected: %s. Released version 3.6 has a breaking bug for haploid mode. More recent versions including nightly builds may work, but are not tested and you must disable preflight checks" % (output))
    except:
    
        # Check for GATK4
        proc = subprocess.Popen(['java', '-jar', jar_path, 'HaplotypeCaller'], stderr=subprocess.PIPE)
        output = proc.stderr.read()
        if "Version:4" in output:
            return 
        else:
            martian.exit("Could not run GATK on %s (%s)\nPossible reasons:\n  1) GATK 3.3 or higher is required\n  2) GATK requires a newer version of Java than the one you have installed" % (hostname, jar_path))


def check_gatk_ref(reference_path):

    if not os.path.exists(os.path.join(reference_path, "fasta", "genome.dict")):
        msg = "GATK requires that you create a .dict reference index file.\n"
        msg += "Note: to use this reference with GATK, also run this Picard command:\n"
        msg += "java -jar /path/to/picard.jar CreateSequenceDictionary R=%s O=%s\n" % (os.path.join(reference_path, "fasta", "genome.fa"), os.path.join(reference_path, "fasta", "genome.dict"))
        msg += "If you have GATK4, use this command:\n"
        msg += "java -jar /path/to/gatk4.jar CreateSequenceDictionary -R=%s" % os.path.join(reference_path, "fasta", "genome.fa")
        martian.exit(msg)
