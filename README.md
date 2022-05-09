# longread_umi_python

## Installation Instructions
The package has not been made executable as of yet. The programs are run directly from `python pipeline.py` at this time.

The program needs to be run in an environment with python version 3.6. Using conda, creating this environment from the command line is accomplished by entering `conda create -n <env name> python=3.6`.

Required packages can be downloaded using the accompanying requirements file. Enter `pip install -r requirements.txt` at the command line once you are in your environment.

Some packages cannot be installed through pip. [This page](https://www.htslib.org/download/) lists three GitHub repositories from which executables are produced. First, run `git clone --recurse-submodules https://github.com/samtools/htslib.git`. Then, for each of them:

1. Clone the github repository (you already have for `htslib`, skip this step for that repository)
2. `cd` into the repository. Enter `make`. Executables should be available in the repository now.
3. Move the executables into the necessary `bin` directory, or change your PATH variable to point to these directories. You will need the `bcftools` and `samtools` executables from their respective directories as well as `bgzip` and `tabix` from the `htslib` directory.

You will also need to download minimap2 via conda. Enter `conda install minimap2`, and if that doesn't work see [here](https://github.com/lh3/minimap2#install) for instructions.
