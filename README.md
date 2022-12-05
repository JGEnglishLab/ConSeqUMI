# longread_umi_python

The package has not been made executable as of yet. The programs are run directly from `python pipeline.py` at this time. 

## Installation Instructions

Note that this program offers three alternative methods for consensus sequence generation, `pairwise`, `lamassemble` and `medaka`. `pairwise` was developed with this package, whereas `lamassemble` and `medaka` were developed elsewhere. A such, the installation instructions will differ depending on which algorithm you plan to use.

Installation instructions assume you have `pip` and `conda` on your computer. If you want all three consensus algorithms available for your environment, follow the Medaka installation instructions first.

### Pairwise Instructions

Required packages can be downloaded using the accompanying requirements file. In your environment enter `pip install -r requirements.txt` to download the necessary dependencies.

### Lamassemble Instructions

Complete 'Pairwise Instructions,' as this method will use the same dependencies.

You will need to download [lamassemble](https://gitlab.com/mcfrith/lamassemble/) separately. Instructions for downloading `lamassemble` can be found in their repository, but for simplicity they are listed here as well. From the command line enter:

`conda install -c bioconda lamassemble`

At present, these instructions download an outdated version of the `lastal` dependency, so you will also need to run:

`conda update last`

Lastly, in accordance with the website instructions, a last-train file needs to be present for `lamassemble` to run properly. From the gitlab repository copy `lamassemble/train/promethion.mat` into this repository's `dependencies_download` folder.

### Medaka Instructions

You will need to download [medaka](https://github.com/nanoporetech/medaka) separately. Instructions for downloading `medaka` can be found in their repository, but for simplicity they are listed here as well. From the command line enter:

`conda create -n medaka -c conda-forge -c bioconda medaka`
(NOTE: This process is significantly faster if you use `mamba` rather than `conda`)

Enter this conda environment by running:

`conda activate medaka`

Complete 'Pairwise Instructions,' as this method will use the same dependencies.

## Finishing Installation

You will need to download [starcode](https://github.com/gui11aume/starcode) using the three steps below.

1. Clone the github repository
2. `cd` into the repository. Enter `make`. Executables should be available in the repository now.
3. Move the executables into the necessary `bin` directory, or change your PATH variable to point to these directories. 
