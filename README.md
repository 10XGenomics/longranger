### Long Ranger OSS release

This repo is the open-source release of Long Ranger 2.2.0. This release is provided for reference and review purposes. It is not yet in a fully buildable state & we are not currently providing support for building and running this code. 
To run Long Ranger, please use the release binaries available at https://support.10xgenomics.com/genome-exome/software/downloads/latest

### Running with docker

Here is an example running [Longranger 2.2.2 at UMCCR](https://hub.docker.com/r/umccr/longranger/) with an Amazon AWS `m5.24xlarge` instance (345 ECUs, 96 vCPUs, 2.5 GHz, Intel Xeon Platinum 8175, 384 GiB memory, EBS only) with an attached 2TB ST1 volume:

```bash
time docker run --privileged -v $PWD/COLO829BL:/data/COLO829BL -v $PWD/refdata-GRCh38-2.1.0:/data/refdata-GRCh38-2.1.0 -v $PWD/pre-called.vcf:/data/pre-called.vcf -v $PWD/output:/data/output -c 95 umccr/longranger:2.2.2 wgs --id="output" --fastqs=/data/COLO829BL --reference=/data/refdata-GRCh38-2.1.0 --jobmode=local --precalled=/data/pre-called.vcf --disable-ui --localcores=95 --localmem=350
```

YMMV, so please adjust the commandline to your machine characteristics. For instance, you might be interested on overriding memory and CPUs per-step with a LongRanger's `overrides.json` file (see **Overriding Default Stage Memory and Thread Requests** in [Using SGE/LSF](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/cluster-mode):

```json
{
"PHASER_SVCALLER_CS.PHASER_SVCALLER._LINKED_READS_ALIGNER.TRIM_READS": { "chunk.threads":24 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._REPORTER.REPORT_COVERAGE": { "join.mem_gb": 12, "chunk.mem_gb": 12, "split.mem_gb": 12 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._REPORTER.REPORT_BASIC": { "split.mem_gb": 12 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._REPORTER.REPORT_SINGLE_PARTITION": { "chunk.mem_gb": 14, "split.mem_gb": 12 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._SNPINDEL_PHASER.ATTACH_PHASING": { "join.mem_gb": 32, "chunk.mem_gb": 32, "split.mem_gb": 32 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._SNPINDEL_PHASER._SNPINDEL_CALLER.POPULATE_INFO_FIELDS": { "chunk.mem_gb": 8, "join.mem_gb": 12 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._SNPINDEL_PHASER._SNPINDEL_CALLER.CALL_SNPINDELS": { "chunk.mem_gb": 12 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._SNPINDEL_PHASER.PHASE_SNPINDELS": { "join.mem_gb": 32 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._LINKED_READS_ALIGNER.MARK_DUPLICATES": { "split.mem_gb": 12 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._LINKED_READS_ALIGNER.MERGE_POS_BAM": { "chunk.mem_gb": 96, "join.mem_gb": 128, "split.mem_gb": 64 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._LINKED_READS_ALIGNER.MERGE_BC_BAM": { "chunk.mem_gb": 96, "join.mem_gb": 128, "split.mem_gb": 64 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._LINKED_READS_ALIGNER._SORT_FASTQ_BY_BARCODE.BUCKET_FASTQ_BY_BC": { "chunk.mem_gb": 12 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._LINKED_READS_ALIGNER._SORT_FASTQ_BY_BARCODE.SORT_FASTQ_BY_BC": { "chunk.mem_gb": 12 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._LINKED_READS_ALIGNER.TRIM_READS": { "chunk.mem_gb": 12 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._LINKED_READS_ALIGNER.BARCODE_AWARE_ALIGNER": { "chunk.mem_gb": 64, "chunk.threads": 32 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._REPORTER.FILTER_BARCODES": { "join.mem_gb": 30 },
"PHASER_SVCALLER_CS.PHASER_SVCALLER._LINKED_READS_ALIGNER._FASTQ_PREP_NEW.BUCKET_FASTQS": { "chunk.threads": 32 }
}
```

In which case the docker commandline will have to include both the `--overrides=/data/overrides.json` flag and `-v $PWD/overrides.json:/data/overrides.json` docker volume mappings.
