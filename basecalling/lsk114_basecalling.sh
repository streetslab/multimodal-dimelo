DORADO_EXE=dorado-1.0.2-linux-x64/bin/dorado
LOCAL_SEQ_DATA=pod5/
EXPERIMENT_ID=experiment_id  # Desired name with which to label all output files
DORADO_MODEL="sup@v5.0.0,5mC_5hmC@v3,6mA@v3"
SEQ_KIT="SQK-LSK114"
BAM_FILE_BASENAME=$EXPERIMENT_ID.dorado_v102_sup_v500_5mC_5hmC_v3_6mA_v3
OUTPUT_DIR=./  # Change if desired
REF_GENOME=chm13v2.0.fasta
FILTER_SCRIPT=dorado_summary_filter.py

$DORADO_EXE basecaller $DORADO_MODEL $LOCAL_SEQ_DATA \
    --modified-bases-threshold 0 \
    --no-trim \
    > $OUTPUT_DIR/$BAM_FILE_BASENAME.bam

$DORADO_EXE trim --sequencing-kit $SEQ_KIT $OUTPUT_DIR/$BAM_FILE_BASENAME.bam > $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.bam
$DORADO_EXE aligner $REF_GENOME $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.bam > $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.bam

$DORADO_EXE summary $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.bam > $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.dorado_summary.tsv
python $FILTER_SCRIPT $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.bam $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.dorado_summary.tsv -t 10

# Optional cleanup
rm $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.bam
