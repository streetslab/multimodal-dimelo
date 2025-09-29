# multimodal-dimelo
multmodal DiMeLo-seq (for paper)

## install

Clone this repo.
```
git clone https://github.com/streetslab/multimodal-dimelo
```

Create your conda environment with the necessary conda dependencies.
```
conda env create -f environment.yml
```

*Note: you can use the --prefix argument to create the environment in a specified directory if you choose. You can then symlink it to your conda home folder.*
```
# OPTIONAL FOR SOME HPC CONFIGURATIONS
conda env create -f environment.yml --prefix /path/to/custom/location
ln -s /path/to/custom/location /home/user/conda/envs/multimodal_dimelo
```

Activate the conda environment.
```
conda activate multimodal_dimelo
```

Install python dependencies and this package with pip. This will bring in the latest dimelo package as well as some other required libraries.
```
pip install -e .
```

Finally, configure this repo to reflect your working environment.
```
cp config.toml.example config.toml
```
Then open `config.toml` and update the empty entries with appropriate paths, etc. for your system. As an example, provide the path to the executable for the UCSC liftOver tool.