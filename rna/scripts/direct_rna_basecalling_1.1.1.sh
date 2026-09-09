#!/bin/bash
#SBATCH --job-name=direct_rna_basecalling_1.1.1
#SBATCH --account=fc_nilah
#SBATCH --partition=savio4_gpu
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:2
#SBATCH --requeue
#SBATCH --time=48:00:00
#SBATCH --output=/clusterfs/nilah/oberon/repos/dimelo-sandbox/direct_rna/basecall_1.1.1_%A.out
#SBATCH --error=/clusterfs/nilah/oberon/repos/dimelo-sandbox/direct_rna/basecall_1.1.1_%A.err

# tar -xzf /global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A.tar.gz -C /global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A
POD5_FOLDER="/global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A/mnt/anna_drive/Projects/AGLemke/RNA004_SampleA"
OUTPUT_BAM="/global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A/basecalled_m6A_pseU.dorado-1.1.1.chm13v2.0.bam"

if [ -f "$OUTPUT_BAM" ]; then
    # Resume from existing file
    /clusterfs/nilah/oberon/software/dorado-1.1.1-linux-x64/dorado-1.1.1-linux-x64/bin/dorado basecaller \
        --reference /clusterfs/nilah/oberon/genomes/chm13v2.0.fa \
        --resume-from "$OUTPUT_BAM" \
        /clusterfs/nilah/oberon/software/dorado-1.1.1-linux-x64/dorado-1.1.1-linux-x64/bin/rna004_130bps_sup@v5.0.0 \
        --modified-bases m6A pseU --estimate-poly-a \
        "$POD5_FOLDER" \
        >> "$OUTPUT_BAM"
else
    # Start fresh
    /clusterfs/nilah/oberon/software/dorado-1.1.1-linux-x64/dorado-1.1.1-linux-x64/bin/dorado basecaller \
        --reference /clusterfs/nilah/oberon/genomes/chm13v2.0.fa \
        /clusterfs/nilah/oberon/software/dorado-1.1.1-linux-x64/dorado-1.1.1-linux-x64/bin/rna004_130bps_sup@v5.0.0 \
        --modified-bases m6A pseU --estimate-poly-a \
        "$POD5_FOLDER" \
        > "$OUTPUT_BAM"
fi

export TMPDIR=/global/scratch/projects/vector_streetslab/oberon/temp
mkdir -p $TMPDIR

source activate whatshap
samtools sort \
    -@ 8 \
    -m 4G \
    -T $TMPDIR/samtools_sort_1.1.1_tmp \
    -o /global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A/basecalled_m6A_pseU_redo.dorado-1.1.1.chm13v2.0.sorted.bam \
    /global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A/basecalled_m6A_pseU_redo.dorado-1.1.1.chm13v2.0.bam
samtools index /global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A/basecalled_m6A_pseU.dorado-1.1.1.chm13v2.0.sorted.bam
