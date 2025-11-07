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
