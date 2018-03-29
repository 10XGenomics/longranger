# Cut out some test events from personalis calls
cut -f1,2,6 personalis_dels.bedpe | grep chr10 > a.bed
awk '{ if ($3 - $2 < 3000) print $0 }' a.bed > test_events.bed
bedtools slop -g hg19.genome -i test_events.bed -b 100 > c.bed
bedtools slop -g hg19.genome -i c.bed -pct -b 0.2 > d.bed
bedtools slop -g hg19.genome -i test_events.bed -b 100 > tp_test_regions.bed
bedtools shuffle -i tp_test_regions.bed -g hg19.genome > fp_test_regions.bed

cat tp_test_regions.bed fp_test_regions.bed > all_regions.bed

# Cut out the corresponding BAM regions to get a modest sized BAM file for testing

rm -f test_regions.bam
rm -f sorted_test_regions.bam

./extract_regions_from_bam.py /mnt/analysis/marsoc/pipestances/HMVT3CCXX/PHASER_SVCALLER_PD/20486/1015.10.0-0/PHASER_SVCALLER_PD/PHASER_SVCALLER/_SNPINDEL_PHASER/ATTACH_PHASING/fork0/files/phased_possorted_bam.bam test_regions.bam all_regions.bed
samtools sort -o sorted_test_regions.bam test_regions.bam 
samtools index sorted_test_regions.bam

rm -f test_regions.bam
rm -f [abcde].bed
rm -f all_regions.bed
