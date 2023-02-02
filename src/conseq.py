import argparse
import os
from Bio import SeqIO
import time

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
        type=OutputDirectory(),
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

def generate_output_name():
    return "ConSeqUMI" + time.strftime("-%Y%m%d-%H%M%S")

class OutputDirectory():

    def __call__(self, name):
        if os.path.isfile(name):
            raise argparse.ArgumentTypeError("The -o or --output argument must be a directory, not a file.")
        if name[-1] != "/": name += "/"
        parentDirectoryPath = "/".join(name.split("/")[:-2])
        if not os.path.isdir(parentDirectoryPath):
            raise argparse.ArgumentTypeError("The -o or --output argument directory requires an existing parent directory.")
        if not os.path.isdir(name):
            name = name[:-1] + "_" + generate_output_name()
        if os.path.isdir(name):
            name += "/" + generate_output_name()
        os.mkdir(name)
        return name
