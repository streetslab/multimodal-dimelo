# Basecalling

Included in this directory are two representative bash scripts for basecalling the ONT sequencing data for this study. One is for single-sample sequencing runs using **SQK-LSK114**, and the other is for multi-sample barcoded sequencing runs using **SQK-NBD114-24**.

Given the raw sequencing files for each experiment, all of the basecalled, modcalled, trimmed, and aligned BAM files used in this study can be exactly recreated using these scripts by adjusting the path variables at the top of the appropriate file, then running the following:

```bash
bash SEQ_KIT_basecalling.sh
```

## Important settings
The following settings are encoded in the scripts themselves and mentioned in the manuscript, but are placed here for ease of reference:
* **dorado version**: `1.0.2`
* **dorado model string**: `sup@v5.0.0,5mC_5hmC@v3,6mA@v3`
* **reference genome**: `chm13v2.0.fasta` (T2T assembly v2.0)
* **filtering settings**:
  * **mean qscore threshold**: `qscore > 10`
  * **minimap2 mapq threshold**: `mapq == 60` (uniquely mapping reads only)

## Direct RNA (m6A / pseU)
The direct-RNA data (RNA004, HEK293T) is processed with Dorado **0.7.2** and samtools by `drna_basecalling.sh`. Before running, download the Dorado 0.7.2 release (linux-x64) from https://github.com/nanoporetech/dorado/releases and unpack it, then download the basecalling model into its `bin/` directory (`cd dorado-0.7.2-linux-x64/bin && ./dorado download --model rna004_130bps_sup@v5.0.0`). The m6A and pseU modification models are fetched automatically next to it on the first run (internet access needed). Edit the path variables at the top of the script (`DORADO_DIR`, `DORADO_MODEL`, `DRNA_DIR`, and the reference genome paths), then run it. It:
* downloads the public POD5 data (~560 GB)
* basecalls with model `rna004_130bps_sup@v5.0.0`, `--modified-bases m6A pseU --estimate-poly-a`, aligned to `chm13v2.0.fa`
* realigns the primary reads (`samtools view -F 0x900`) with `--mm2-preset splice:hq` to the hs1 + rDNA (`chrR`) genome from George, Pimkin & Paralkar, *JBC* 2023. Clone https://github.com/vikramparalkar/rDNA-Mapping-Genomes and use `Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa` (extract the `.tar.gz` if needed).

The final BAM (`basecalled_m6A_pseU.splice.hs1-rDNA_v1.0.sorted.bam`) is the one used by `direct_rna.ipynb`.
