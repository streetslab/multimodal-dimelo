# Direct RNA (RNA004) basecalling with m6A + pseU modification calling, on the public HEK293T "Sample A" data
# (Hewel et al., https://www.biorxiv.org/content/10.1101/2024.07.25.605188v3.full.pdf)
# Stages: download, basecall (aligned to chm13), then spliced realignment of the primary reads to hs1 + rDNA.
# Single-pass alternative: pass --reference $RDNA_GENOME --mm2-preset splice:hq to the basecaller and skip the realignment stage;
# this should give essentially identical results but was not the pipeline used for the manuscript.
DORADO_DIR=dorado-0.7.2-linux-x64
DORADO_EXE=$DORADO_DIR/bin/dorado
DORADO_MODEL=$DORADO_DIR/bin/rna004_130bps_sup@v5.0.0  # model directory as downloaded alongside dorado
DRNA_URL=ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR152/ERR15278635/RNA004_HEK293T_A.tar.gz  # POD5 tarball is large (~560 GB)
DRNA_DIR=dRNA/RNA004_HEK293T_A  # Change if desired
LOCAL_SEQ_DATA=$DRNA_DIR/mnt/anna_drive/Projects/AGLemke/RNA004_SampleA
REF_GENOME=chm13v2.0.fa  # intermediate alignment target only; the final alignment is to RDNA_GENOME below
# hs1 (T2T-CHM13) + chrR rDNA contig (https://github.com/vikramparalkar/rDNA-Mapping-Genomes)
RDNA_GENOME=hs1-rDNA_v1.0.fa
BAM_FILE_BASENAME=$DRNA_DIR/basecalled_m6A_pseU_redo.chm13v2.0
OUTPUT_BASENAME=$DRNA_DIR/basecalled_m6A_pseU.splice.hs1-rDNA_v1.0
TMPDIR=./tmp  # Change if desired
mkdir -p $TMPDIR

if [ ! -d $LOCAL_SEQ_DATA ]; then
    mkdir -p $DRNA_DIR
    wget -P $DRNA_DIR $DRNA_URL
    tar -xzf $DRNA_DIR/RNA004_HEK293T_A.tar.gz -C $DRNA_DIR
fi

$DORADO_EXE basecaller \
    --reference $REF_GENOME \
    $DORADO_MODEL \
    --modified-bases m6A pseU --estimate-poly-a \
    $LOCAL_SEQ_DATA \
    > $BAM_FILE_BASENAME.bam

samtools sort -@ 8 -m 4G -T $TMPDIR/samtools_sort_tmp -o $BAM_FILE_BASENAME.sorted.bam $BAM_FILE_BASENAME.bam
samtools index $BAM_FILE_BASENAME.sorted.bam

# Spliced realignment to hs1 + rDNA; only primary records are realigned (-F 0x900) so each read appears once.
# The final BAM ($OUTPUT_BASENAME.sorted.bam) is the one read by direct_rna.ipynb
samtools view -b -F 0x900 -@ 8 -o $OUTPUT_BASENAME.primary_input.bam $BAM_FILE_BASENAME.sorted.bam

$DORADO_EXE aligner $RDNA_GENOME $OUTPUT_BASENAME.primary_input.bam --mm2-preset splice:hq > $OUTPUT_BASENAME.bam

samtools sort -@ 8 -m 4G -T $TMPDIR/samtools_sort_realign -o $OUTPUT_BASENAME.sorted.bam $OUTPUT_BASENAME.bam
samtools index $OUTPUT_BASENAME.sorted.bam
