# Direct RNA / ribosomal RNA base modification figures

Handoff copy of the scripts and notebook used to generate the paper's direct-RNA
(m6A / pseudouridine) figures over ribosomal RNA (18S/28S/5.8S rRNA) and other
transcripts. Paths inside the scripts/notebook are left exactly as they were run
on Savio — none of this has been rewired into the reproducibility repo's path
management yet, so file paths still point at `/clusterfs/nilah/oberon/...` and
`/global/scratch/projects/vector_streetslab/oberon/...`. Originals live at
`/clusterfs/nilah/oberon/repos/dimelo-sandbox/direct_rna/`.

## 1. Where the data came from

Raw direct-RNA nanopore data (RNA004 kit, HEK293T cells, "Sample A") was
downloaded from ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR152/ERR15278635/RNA004_HEK293T_A.tar.gz as a
tarball (per https://www.biorxiv.org/content/10.1101/2024.07.25.605188v3.full.pdf) and extracted on Savio scratch:

```
/global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A.tar.gz
  -> extracted to ->
/global/scratch/projects/vector_streetslab/oberon/dRNA/RNA004_HEK293T_A/mnt/anna_drive/Projects/AGLemke/RNA004_SampleA/
```

That folder contains the POD5 signal files that basecalling reads directly.
This raw data is **not** copied here (it's ~large POD5 data on scratch) — only
the processing code is included. 

## 2. Processing pipeline

1. **`scripts/direct_rna_basecalling_1.1.1.sh`** — Dorado basecalling of the
   POD5 folder above.
   - Model: `rna004_130bps_sup@v5.0.0`
   - Modified bases called: `m6A`, `pseU` (`--modified-bases m6A pseU`), plus
     `--estimate-poly-a`
   - Aligned during basecalling to T2T-CHM13 (`/clusterfs/nilah/oberon/genomes/chm13v2.0.fa`)
   - Two Dorado versions were tried (0.7.2 and 1.1.1) — **1.1.1 is the one
     actually carried forward** into the realignment step and the figures
     (see `direct_rna_realign_to_rDNA.sh` below, which reads the `chm13v2.0`
     sorted BAM produced by the 0.7.2/1.1.1 script family under the
     `basecalled_m6A_pseU_redo.chm13v2.0.sorted.bam` name).
   - Dorado binaries used:
     `/clusterfs/nilah/oberon/software/dorado-0.7.2-linux-x64/`
     `/clusterfs/nilah/oberon/software/dorado-1.1.1-linux-x64/`
   - Output sorted/indexed with `samtools` (see environment below).

2. **`scripts/direct_rna_realign_to_rDNA.sh`** — realigns the CHM13-mapped BAM
   with `dorado aligner` (`--mm2-opts "-x splice:hq"`, i.e. spliced long-read
   alignment) against a **custom ribosomal DNA reference genome**:
   `/clusterfs/nilah/oberon/repos/rDNA-Mapping-Genomes/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa`.
   - This reference is the hs1 (T2T-CHM13) genome with an added `chrR`
     ribosomal-DNA-repeat contig, from George, Pimkin & Paralkar,
     "Construction and validation of customized genomes for human and mouse
     ribosomal DNA mapping," *JBC* (2023):
     https://www.jbc.org/article/S0021-9258(23)01794-5/fulltext.
     Publicly deposited copy used here:
     `/clusterfs/nilah/oberon/repos/rDNA-Mapping-Genomes` (separate git repo,
     includes `.bed` annotations for rRNA/IGS regions on `chrR`).
   - Output: `basecalled_m6A_pseU.splice.hs1-rDNA_v1.0.sorted.bam` — this is
     the BAM the notebook actually reads.

3. **`scripts/stop_codon_from_annotation.py`** — a standalone Python script
   (no Slurm wrapper; run directly with the `whatshap` env's Python) that
   parses `hs1.ncbiRefSeq.gtf.gz` for high-confidence (`NM_`/`NR_`) curated
   transcripts and writes stop-codon positions to a BED file:
   `/clusterfs/nilah/oberon/datasets/annotations/hs1_stop_codons_high_confidence.bed`.
   This BED is one of the region inputs to the stop-codon enrichment-profile
   figure in the notebook (`m6A_pseU_profile_around_stop_codons.svg`).

4. **`notebooks/direct_rna.ipynb`** — the actual figure-generation notebook.
   Uses `dimelo-toolkit`'s `parse_bam.pileup` to build a modification pileup
   from the realigned BAM against the `hs1-rDNA_v1.0.fa` reference, then plots
   with `plot_enrichment`, `plot_enrichment_profile`, and `plot_read_browser`.
   Key regions plotted (all on the custom `chrR` contig):
   - `RRNA_18S = "chrR:12996-14864"`
   - `RRNA_28S = "chrR:17259-22309"`
   - `RRNA_5_8S = "chrR:15935-16091"`
   - `RRNA_all = "chrR:0-44837"`
   Modification motifs: `A,0,a` (m6A) and `T,0,17802` (pseudouridine, Ψ), plus
   sequence-context variants (`DRACH,2,a`, `GTTC,2,17802`, `TSTAG,2,17802`, …).

   Note: the notebook as saved has several early exploratory cells commented
   out (alternate pileups/regions that were tried and abandoned); the cells
   that actually produced the saved figures are the ones left uncommented,
   in particular the `chr1_R_pileup` pileup call and the `plot_enrichment` /
   `plot_enrichment_profile` / `plot_read_browser` calls further down.

## 3. Figures

Copied as-generated into `figures/`:

| File | Content |
|---|---|
| `18SRNA_flanking_read_browser.svg` | Read browser, 18S rRNA + flanking region |
| `18SRNA_spanning_read_browser.svg` | Read browser, reads spanning full 18S rRNA |
| `m6A_pseU_enrichment_across_18S ribosomal subunit.svg` | m6A/Ψ enrichment, 18S subunit |
| `m6A_pseU_enrichment_across_18S_rRNA_gene_body.svg` | m6A/Ψ enrichment, 18S gene body |
| `m6A_pseU_enrichment_across_ribosomal_subunits.svg` | m6A/Ψ enrichment, all rRNA subunits |
| `m6A_pseU_enrichment_across_chr1_transcripts.svg` | m6A/Ψ enrichment, chr1 transcripts (comparison/control) |
| `m6A_pseU_profile_around_stop_codons.svg` | m6A/Ψ metagene profile around annotated stop codons |

## 4. Environment setup

Two separate conda environments were used for two separate stages:

- **Basecalling / alignment / samtools sort+index** — conda env `whatshap`
  (`/clusterfs/nilah/oberon/environments/whatshap`, samtools 1.21 / htslib
  1.21). Activated in the Slurm scripts via `source activate whatshap`.
  Dorado itself is a standalone binary, not a conda package:
  - Dorado 0.7.2: `/clusterfs/nilah/oberon/software/dorado-0.7.2-linux-x64/`
  - Dorado 1.1.1: `/clusterfs/nilah/oberon/software/dorado-1.1.1-linux-x64/`
    (reports version `1.1.1+e72f1492`)

- **Notebook / plotting** — conda env `dimelo-toolkit`

## 5. Supporting reference data (not copied here)

- `/clusterfs/nilah/oberon/repos/rDNA-Mapping-Genomes/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa`
  — custom rDNA-mapping reference genome (separate git repo, see its own
  README for citation/provenance).
- `/clusterfs/nilah/oberon/genomes/chm13v2.0.fa` — standard T2T-CHM13v2.0,
  used only as the intermediate basecalling-time alignment target.
- `/clusterfs/nilah/oberon/datasets/annotations/hs1.ncbiRefSeq.gtf.gz` — RefSeq
  annotation for hs1, input to `stop_codon_from_annotation.py`.
- `/clusterfs/nilah/oberon/datasets/annotations/hs1_stop_codons_high_confidence.bed`
  — generated by `stop_codon_from_annotation.py`, not regenerated automatically
  by anything else in this folder.
