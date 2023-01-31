import argparse

def set_command_line_settings():
    parser = argparse.ArgumentParser(description="")
    commandParser = parser.add_subparsers(dest="command", help="ConSeq Functions")
    commandParser.add_parser(
        "umi",
        help="Extracts UMIs and target sequences from .fastq files",
    )
    return parser

def check_for_invalid_input(parser, args): pass
