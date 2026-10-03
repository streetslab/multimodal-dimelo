OUTPUT_DIR=./  # Change if desired
REF_GENOME=chm13v2.0.fa
VCF_FILE=NA12878.t2t_chm13v2.vcf.gz  # Change to reflect location of the VCF file
INPUT_BAM=path/to/file_name.bam  # Change to reflect target file
# AUTO POPULATED VARIABLES
UNTAGGED_BAM=$OUTPUT_DIR/untagged.bam
BAM_ID=$(basename $INPUT_BAM .bam)
HAPLOTAGGED_BAM=$OUTPUT_DIR/$BAM_ID.haplotagged.bam

whatshap haplotag --ignore-read-groups --reference $REF_GENOME $VCF_FILE $INPUT_BAM | samtools view -@16 -b -o $HAPLOTAGGED_BAM
samtools split -u $UNTAGGED_BAM -d HP -@16 $HAPLOTAGGED_BAM
samtools index -@20 $BAM_ID.haplotagged_1.bam
samtools index -@20 $BAM_ID.haplotagged_2.bam
mv $BAM_ID.haplotagged_1.bam* $OUTPUT_DIR
mv $BAM_ID.haplotagged_2.bam* $OUTPUT_DIR