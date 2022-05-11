# DepMap autoencoder normalization

This repository contains R and Python scripts to train autoencoders on a processed version of the 20Q2 DepMap dataset. 

To better train any sort of deep learning algorithm, the DepMap was processed before training in the following order:
- Row-standardization
- Clipping
- Min-max scaling

## Setup

### Dependencies and versioning

Python package dependencies:
- PyTorch (1.6.0)
- NumPy (1.17.4)
- Pandas (0.25.1)

R package dependencies:
- ggplot2 (3.3.3)
- ggthemes (4.2.4)
- pheatmap (1.0.12)
- argparse (2.0.3)
- gplots (3.1.1)
- stringi (1.5.3)
- devtools (2.3.2)
- RColorBrewer (1.1.2)
- FLEX (1.0.0, available at https://github.com/csbio/FLEX_R)

All R code tested only with version 3.6.3 on Ubuntu 20.04.3 LTS. All Python code tested only with Python version 3.7.3 on Ubuntu 20.04.3 LTS. There may or may not be compatability issues with packages built on R version 4.0.0 and later. 

### Data setup

The original, processed 20Q2 DepMap dataset must be placed within the input subdirectory. This is downloadable <here>. 

Alternately, a given release of the DepMap may be downloadable from the DepMap website here: https://depmap.org/portal/download/all/. This may then be processed with the script <scriptname.R> and placed in the input subdirectory.

## Workflow

To run the pipeline, navigate to this directory and call the train_autoencoder.sh bash script: `./train_autoencoder.sh`