# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import tenkit.log_subprocess
import tenkit.bio_io as tk_io
import tenkit.tabix as tk_tabix

def canonicalize(filename, output_name):
        with open(output_name+"tmp.vcf",'w') as canon_file:
            # Left-align indels and output variants as constitutent indels
            tenkit.log_subprocess.check_call(['vcfallelicprimitives','--keep-info', '-t','VCFALLELICPRIMITIVE', filename], stdout=canon_file)
        with open(output_name+"tmp2.vcf",'w') as fixed_vcf:
            # the reason we are doing this is because vcfallelicprimitives screws up the vcf format in some places in the info fields
            #  it changes some tags from KEY=.; to KEY=; which is invalid. bcftools fixes this, but we dont actually want to filter anything
            #   this should do that
            tenkit.log_subprocess.check_call(['bcftools', 'filter', output_name+"tmp.vcf"],stdout=fixed_vcf)
        with open(output_name+"tmp3.vcf",'w') as unphased_file:
            vcf_in = tk_io.VariantFileReader(output_name+"tmp2.vcf")
            unsupported_genotype_filter = [("UNSUPPORTED_GENOTYPE", "If genotype field contains '.' we assume that this is due to making a single sample vcf from a multiple sample vcf in which this sample does not contain the variant.")]
            tenx = ('TENX','0','Flag',"called by 10X", None, None)
            vcf_out = tk_io.VariantFileWriter(unphased_file, template_file=open(output_name+"tmp2.vcf"), new_info_fields=[tenx], new_filters = unsupported_genotype_filter)
            for record in vcf_in.record_getter():
                sample = record.samples[0]
                unsupported_genotype = False
                try:
                    if len(sample.gt_alleles) > 0:
                        genotype1 = sample.gt_alleles[0]
                        if genotype1 == '.':
                            unsupported_genotype = True
                    else:
                        unsupported_genotype = True
                    if len(sample.gt_alleles) > 1:
                        genotype2 = sample.gt_alleles[1]
                        if genotype2 is '.':
                            unsupported_genotype = True
                except:
                    unsupported_genotype = True
                if unsupported_genotype:
                    record.FILTER = ["UNSUPPORTED_GENOTYPE"]
                vcf_out.write_record(record)
        tk_tabix.sort_vcf(output_name+"tmp3.vcf", output_name)
