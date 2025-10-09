# multimodal-dimelo
multmodal DiMeLo-seq (for paper)

## Installation

Clone this repo:
```
git clone https://github.com/streetslab/multimodal-dimelo
```

Create your conda environment with the necessary conda dependencies:
```
conda env create -f environment.yml
```
*Note: Instead of installing the environment in this directory, you can use the --prefix argument to create the environment in a specified location. You can then symlink it to your conda home folder. This may be useful for installation on cluster systems with limited user permissions.*
```
# OPTIONAL FOR SOME HPC CONFIGURATIONS
conda env create -f environment.yml --prefix /PATH/TO/CUSTOM/LOCATION
ln -s /PATH/TO/CUSTOM/LOCATION /home/user/conda/envs/multimodal_dimelo
```

Activate the conda environment:
```
conda activate multimodal_dimelo
```

Install python dependencies and this package with pip. This will bring in the latest dimelo package as well as some other required libraries:
```
pip install .
```

### Installation and configuration of external tools
Download and/or install the following:
* [UCSC liftOver tool](https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [samtools](https://www.htslib.org/)


To configure this repository to reflect your working environment, first copy the example configuration file:
```
cp config.toml.example config.toml
```
Then open `config.toml` and update the empty entries with appropriate paths to the executables for your system (e.g. `/usr/bin/bedtools`)

# Running analysis

TODO: This needs to include basecalling and/or data download instructions for the BAMs, maybe?

To regenerate the figures for this manuscript, first run third_party_data.ipynb to download all necessary data and perform necessary processing to generate reference files (e.g. BEDs).

Then run any of the analysis notebooks in any order:
1. CTCF_analysis.ipynb
2. LMNB1_analysis.ipynb