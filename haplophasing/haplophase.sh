OUTPUT_DIR=./  # Change if desired
REF_GENOME=chm13v2.0.fasta
VCF_FILE=NA12878.chm13v2.vcf.gz
INPUT_BAM=path/to/file_name.bam  # Change to reflect target file
OUTPUT_ID=output_id  # Change to a useful ID for the output bam files

whatshap haplotag --ignore-read-groups --output-haplotag-list $OUTPUT_DIR/haplotagged.tsv --reference $REF_GENOME $VCF_FILE $INPUT_BAM | samtools view -@8 -b -o haplotagged.bam
whatshap split --output-h1 $OUTPUT_DIR/$OUTPUT_ID.h1.bam --output-h2 $OUTPUT_DIR/$OUTPUT_ID.h2.new.bam $OUTPUT_DIR/haplotagged.bam $OUTPUT_DIR/haplotagged.tsv
rm $OUTPUT_DIR/haplotagged.bam
rm $OUTPUT_DIR/haplotagged.tsv
samtools index -@20 $OUTPUT_DIR/$OUTPUT_ID.h1.bam
samtools index -@20 $OUTPUT_DIR/$OUTPUT_ID.h2.bam