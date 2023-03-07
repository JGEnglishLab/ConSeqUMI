# ConseqUMI

## Basic Installation

### Conda environment

It is recommended you create an [anaconda](https://docs.anaconda.com/anaconda/install/) environment with python for installing all the necessary packages. Presuming you have `anaconda` installed, you can create a python environment with the following commands from the command line:

```
conda create -n conseq
conda activate conseq
conda install python
```

### Package dependencies

Clone this repository and download the dependencies with the following commands.

```
git clone https://github.com/JGEnglishLab/ConseqUMI.git
cd ConseqUMI
pip install .
```

## Starcode Installation

You will need to download [starcode](https://github.com/gui11aume/starcode) using the three steps below.

1. Clone the starcode github repository.
2. `cd` into the repository. Enter `make`. Executables should be available in the repository now.
3. Move the executables into the necessary `bin` directory, or change your PATH variable to point to these directories. 

## Consensus Sequence Instructions

Note that this program offers three alternative methods for consensus sequence generation, `pairwise`, `lamassemble` and `medaka`. `pairwise` is a custom consensus method developed with this package, whereas `lamassemble` and `medaka` were developed elsewhere. A such, the installation instructions will differ depending on which algorithm you plan to use.

Installation instructions assume you have `pip` and `conda` on your computer. If you want all three consensus algorithms available for your environment, follow the Medaka installation instructions first.

### Pairwise Instructions

No extra instructions are necessary for the `pairwise` method.

### Lamassemble Instructions

Complete 'Pairwise Instructions,' as this method will use the same dependencies.

You will need to download [lamassemble](https://gitlab.com/mcfrith/lamassemble/) separately. Instructions for downloading `lamassemble` can be found in their repository, but for simplicity they are listed here as well. From the command line enter:

`conda install -c bioconda lamassemble`

At present, these instructions download an outdated version of the `lastal` dependency, so you will also need to run:

`conda update last`

Lastly, in accordance with the website instructions, a last-train file needs to be present for `lamassemble` to run properly. From the gitlab repository copy `lamassemble/train/promethion.mat` into this repository's `dependencies_download` folder.

To adjust the default lamassemble settings used by this program, edit the `lamassembleCommandLine` variable in the `consensus_generators/config.py` file.

### Medaka Instructions

You will need to download [medaka](https://github.com/nanoporetech/medaka) separately. Instructions for downloading `medaka` can be found in their repository, but for simplicity they are listed here as well. From the command line enter:

`conda create -n medaka -c conda-forge -c bioconda medaka`
(NOTE: This process is significantly faster if you use `mamba` rather than `conda`)

Enter this conda environment by running:

`conda activate medaka`

Complete 'Pairwise Instructions,' as this method will use the same dependencies.

To adjust the default medaka settings used by this program, edit the `medakaCommandLine` variable in the `consensus_generators/config.py` file.

## Running ConseqUMI

Run `conseq` to start a user-friendly GUI. For command line instructions, run `conseq -h`.
