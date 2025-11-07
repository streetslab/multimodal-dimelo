DORADO_EXE=dorado-1.0.2-linux-x64/bin/dorado
LOCAL_SEQ_DATA=pod5/
MINKNOW_ID=c65cb8c5-4595-4682-98e1-90cd57cbeded  # Unique ID given to each sequencing run by MinKNOW; replace with correct ID
EXPERIMENT_ID=experiment_id  # Desired name with which to label all output files
DORADO_MODEL="sup@v5.0.0,5mC_5hmC@v3,6mA@v3"
SEQ_KIT="SQK-NBD114-24"
BAM_FILE_BASENAME=$EXPERIMENT_ID.dorado_v102_sup_v500_5mC_5hmC_v3_6mA_v3
OUTPUT_DIR=./  # Change if desired
REF_GENOME=chm13v2.0.fasta
EXPERIMENT_BARCODES=(01 02 03 04)  # Update to accurately reflect barcodes used in this experiment
FILTER_SCRIPT=dorado_summary_filter.py

$DORADO_EXE basecaller $DORADO_MODEL $LOCAL_SEQ_DATA \
    --kit-name $SEQ_KIT \
    --modified-bases-threshold 0 \
    --no-trim \
    > $OUTPUT_DIR/$BAM_FILE_BASENAME.bam

$DORADO_EXE trim --sequencing-kit $SEQ_KIT $OUTPUT_DIR/$BAM_FILE_BASENAME.bam > $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.bam
$DORADO_EXE aligner $REF_GENOME $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.bam > $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.bam

$DORADO_EXE summary $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.bam > $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.dorado_summary.tsv
python $FILTER_SCRIPT $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.bam $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.dorado_summary.tsv -t 10

$DORADO_EXE demux --output-dir $OUTPUT_DIR/demuxed --no-classify $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.filtered.q10.bam

for barcode in ${EXPERIMENT_BARCODES[@]}; do
    demuxed_bam=$OUTPUT_DIR/demuxed/${MINKNOW_ID}_${SEQ_KIT}_barcode$barcode.bam
    sorted_bam=$OUTPUT_DIR/demuxed/barcode$barcode.sorted.bam
    # samtools sort and index
    samtools sort -m1g -@20 -o $sorted_bam $demuxed_bam
    samtools index -@20 $sorted_bam
done

# Optional cleanup
rm $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.bam
rm $OUTPUT_DIR/$BAM_FILE_BASENAME.trim.align.filtered.q10.bam
rm $OUTPUT_DIR/demuxed/${MINKNOW_ID}_*.bam
