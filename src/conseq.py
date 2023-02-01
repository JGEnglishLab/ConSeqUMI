import argparse
import os
from Bio import SeqIO

def set_command_line_settings():
    parser = argparse.ArgumentParser(description="")
    commandParser = parser.add_subparsers(dest="command", help="ConSeq Functions")
    umiParser = commandParser.add_parser(
        "umi",
        help="Extracts UMIs and target sequences from .fastq files",
    )
    umiParser.add_argument(
        "-i",
        "--input",
        type=InputDirectory(),
        required=True,
        help="Path to folder that only contains input Nanopore read fastq files.",
    )
    umiParser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path for folder output. Folder should not currently exist.",
    )
    umiParser.add_argument(
        "-a",
        "--adapters",
        type=str,
        required=True,
        help="A text file with f, F, r, R adapters listed. Defaults to: GAGTGTGGCTCTTCGGAT, ATCTCTACGGTGGTCCTAAATAGT, AATGATACGGCGACCACCGAGATC, and CGACATCGAGGTGCCAAAC, respectively.",
    )

    return parser

class InputDirectory():
    def __init__(self):
        self.allowedFileTypes = set(["fastq","fq"])
    def __call__(self, name):
        if os.path.isfile(name):
            raise argparse.ArgumentTypeError("The -i or --input argument must be a directory, not a file.")
        if name[-1] != "/": name += "/"
        if not os.path.isdir(name):
            raise argparse.ArgumentTypeError("The -i or --input argument must be an existing directory.")
        files = os.listdir(name)
        if len(files) == 0:
            raise argparse.ArgumentTypeError("The -i or --input argument directory must not be empty.")
        records = []
        for file in files:
            if file.split('.')[-1] not in self.allowedFileTypes:
                raise argparse.ArgumentTypeError("The -i or --input argument directory must only contain fastq files (.fq or .fastq).")
            else:
                records.extend([record for record in SeqIO.parse(name + file, "fastq")])
        return records
