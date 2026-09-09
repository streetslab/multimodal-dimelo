#!/bin/bash
#SBATCH --job-name=direct_rna_realign_to_rDNA
#SBATCH --account=fc_nilah
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --output=/clusterfs/nilah/oberon/repos/dimelo-sandbox/direct_rna/slurm_outputs/realign_rDNA_%A.out
#SBATCH --error=/clusterfs/nilah/oberon/repos/dimelo-sandbox/direct_rna/slurm_outputs/realign_rDNA_%A.err

/clusterfs/nilah/oberon/software/dorado-1.1.1-linux-x64/dorado-1.1.1-linux-x64/bin/dorado aligner \
/clusterfs/nilah/oberon/repos/rDNA-Mapping-Genomes/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa \
/global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A/basecalled_m6A_pseU_redo.chm13v2.0.sorted.bam \
--mm2-opts "-x splice:hq" \
> /global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A/basecalled_m6A_pseU.splice.hs1-rDNA_v1.0.bam

export TMPDIR=/global/scratch/projects/vector_streetslab/oberon/temp
mkdir -p $TMPDIR

source activate whatshap
samtools sort \
    -@ 8 \
    -m 4G \
    -T $TMPDIR/samtools_sort_realign \
    -o /global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A/basecalled_m6A_pseU.splice.hs1-rDNA_v1.0.sorted.bam \
    /global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A/basecalled_m6A_pseU.splice.hs1-rDNA_v1.0.bam
samtools index /global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A/basecalled_m6A_pseU.splice.hs1-rDNA_v1.0.sorted.bam
