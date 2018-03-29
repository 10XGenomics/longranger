# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import os
import tenkit.log_subprocess
import subprocess
from subprocess import CalledProcessError


def check_gatk_version(gatk_path):

    try:
        output = subprocess.check_output(["java", "-jar", gatk_path, "--version"])
        is_gatk3 = "3.2" in output or "3.3" in output or "3.4" in output or "3.5" in output or "3.7" in output or "3.8" in output or "3.9" in output
        is_gatk37 = "3.7" in output or "3.8" in output or "3.9" in output
        if is_gatk3:
            return ("3", is_gatk37)

    except CalledProcessError: 
        pass


    proc = subprocess.Popen(['java', '-jar', gatk_path, 'HaplotypeCaller'], stderr=subprocess.PIPE)
    output = proc.stderr.read()
    if "Version:4" in output:
        return ("4", True)

    return None

def run_variant_caller(variant_caller, gatk_path, mem_gb, fasta_path, bamfile, output, bedfile=None, haploid_mode=False):


    if variant_caller == 'freebayes':

        with open(output, 'w') as output_file:
            args = ['freebayes', '-f', fasta_path, '-b', bamfile, '-0']
            if bedfile is not None:
                args.extend(['-t', bedfile])

            if haploid_mode:
                args.extend(['-p', '1'])

            print args
            tenkit.log_subprocess.check_call(args, stdout=output_file)


    elif variant_caller == 'gatk':

        gatk_version = check_gatk_version(gatk_path)
        if gatk_version is None:
            raise RuntimeError("Couldn't determine GATK version -- check you GATK path")
    
        (main_ver, is_gatk37) = gatk_version

        jre_args = ["java", "-Xmx%dG" % int(mem_gb), '-XX:ParallelGCThreads=2']
    
        if os.environ.get("TMPDIR"):
            jre_args.append("-Djava.io.tmpdir=%s" % os.environ["TMPDIR"])

        base_args = jre_args + ['-jar', gatk_path]
        if main_ver == "3":
            base_args.extend(["-T", "HaplotypeCaller"])
        elif main_ver == "4":
            base_args.extend(["HaplotypeCaller"])
        

        if is_gatk37:
            emit_args = []
        else:
            emit_args = ['-stand_emit_conf', '10']

        if main_ver == "3":
            mode_args = ['--genotyping_mode', 'DISCOVERY', '-stand_call_conf', '30']
        else:
            mode_args = ['--genotyping-mode', 'DISCOVERY']

        main_args = ['-R', fasta_path, '-I', bamfile] + mode_args + emit_args 

        if bedfile is not None:
            main_args += ['-L', bedfile]

        if main_ver == "3":
            main_args += ['-o', output]
        else:
            main_args += ['-O', output] 

        args = base_args + main_args
        print args
        tenkit.log_subprocess.check_call(args)

    else:
        raise RuntimeError('Variant caller not supported: ' + variant_caller)
