import argparse

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
        type=str,
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

def check_for_invalid_input(parser, args): pass
