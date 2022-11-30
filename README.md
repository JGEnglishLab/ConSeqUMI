# longread_umi_python

## Installation Instructions

Required packages can be downloaded using the accompanying requirements file. Enter `pip install -r requirements.txt` at the command line once you are in your environment. If you plan on using the `medaka` consensus option, you will need to load different dependencies. Enter `pip install -r requirements_medaka.txt` instead.

If you plan on using the `medaka` or `lamassemble` consensus options, you will need to download these respective packages separately.
[medaka link](https://github.com/nanoporetech/medaka)
[lamassemble link](https://gitlab.com/mcfrith/lamassemble/)

The package has not been made executable as of yet. The programs are run directly from `python pipeline.py` at this time.

### Finishing Installation

Some packages cannot be installed through pip or conda. [This page](https://www.htslib.org/download/) lists three GitHub repositories from which executables are produced. First, run `git clone --recurse-submodules https://github.com/samtools/htslib.git`. Then, for each of them:

1. Clone the github repository
2. `cd` into the repository. Enter `make`. Executables should be available in the repository now.
3. Move the executables into the necessary `bin` directory, or change your PATH variable to point to these directories. 

You will also need to download [starcode](https://github.com/gui11aume/starcode) using the three steps above.
