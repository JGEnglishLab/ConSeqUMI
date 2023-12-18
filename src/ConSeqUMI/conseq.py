import argparse
import os
from Bio import SeqIO
import time
from ConSeqUMI.umi import umi
from ConSeqUMI.gui import gui
from ConSeqUMI.consensus import benchmark, consensus
from ConSeqUMI.Printer import Printer
from shutil import which


def main():
    printer = Printer()
    printer("parsing and loading arguments")
    parser = set_command_line_settings()
    args = vars(parser.parse_args())
    if args["command"] == "gui" or not args["command"]:
        gui.main()
    printer(f"output directory: {args['output']}")
    if args["command"] == "umi":
        umi.main(args)
    if args["command"] == "cons":
        consensus.main(args)
    if args["command"] == "benchmark":
        benchmark.main(args)



def set_command_line_settings():
    parser = argparse.ArgumentParser(description="")
    commandParser = parser.add_subparsers(dest="command", help="ConSeq Functions")
    guiParser = commandParser.add_parser(
        "gui",
        help="Starts the program GUI.",
    )
    umiParser = commandParser.add_parser(
        "umi",
        help="Extracts UMIs and target sequences from .fastq files",
    )
    umiParser.add_argument(
        "-i",
        "--input",
        type=InputDirectory("umi"),
        required=True,
        help="Path to folder that only contains input Nanopore read fastq files.",
    )
    umiParser.add_argument(
        "-o",
        "--output",
        type=OutputDirectory("umi"),
        required=True,
        help="Path for folder output. Folder should not currently exist.",
    )
    umiParser.add_argument(
        "-a",
        "--adapters",
        type=AdapterFile(),
        required=True,
        help="A text file with f, F, r, R adapters listed in order.",
    )
    umiParser.add_argument(
        "-u",
        "--umiLength",
        type=ConseqInt("umiLength"),
        default=0,
        help="The expected length of any UMI found, minimum 10. Providing this option loosens front adapter requirements and generally results in increased quantity of UMIs and target sequences found.",
    )
    consParser = commandParser.add_parser(
        "cons",
        help="Finds a consensus sequence for each fastq file in a given directory and writes them to a single output fasta file.",
    )
    consParser.add_argument(
        "-i",
        "--input",
        type=InputDirectory("cons"),
        required=True,
        help="Path to directory that only contains fastq files. Note that each individual fastq file should contain sequences that contribute to a single consensus. If directing at the 'umi' command output, this will be the 'bins' directory in the 'umi' command output.",
    )
    consParser.add_argument(
        "-o",
        "--output",
        type=OutputDirectory("consensus"),
        required=True,
        help="Path for folder output. The folder will be created with a time stamp.",
    )
    consParser.add_argument(
        "-c",
        "--consensusAlgorithm",
        type=ConsensusAlgorithmText(),
        default="pairwise",
        help="An option between two consensus sequence algorithms. Default is a customized algorithm that relies on pairwise alignment, which can be slow for larger sequences. Options: pairwise (default), lamassemble",
    )
    consParser.add_argument(
        "-m",
        "--minimumReads",
        type=ConseqInt("minimumReads"),
        default=50,
        help="Minimum number of cluster reads required to generate a consensus sequence. Default is 50.",
    )
    consParser.add_argument(
        "-p",
        "--processNum",
        type=ConseqInt("processNum"),
        default=1,
        help="Number of processes to run. By default it will only use 1. If you enter a number beyond the number of processes your computer is capable of, the number of processes will automatically be set to the maximum level for your computer.",
    )
    consParser.add_argument(
        "-l",
        "--lastTrain",
        type=LastTrainFile(),
        help="Path to a last-train mat file for lamassemble. If you have already put the path to the desired file in the consensus/config.py file, this flag is unnecessary and should not be used.",
    )
    benchmarkParser = commandParser.add_parser(
        "benchmark",
        help="Creates a benchmarking data analysis file for evaluating the accuracy of a provided consensus sequence algorithm when applied to a given input fastq file.",
    )
    benchmarkParser.add_argument(
        "-i",
        "--input",
        type=InputFile("input"),
        required=True,
        help="Path to a fastq file. Note that the fastq file should contain sequences that contribute to a single consensus. If directing at the 'umi' command output, this will be in the 'bins' directory in the 'umi' command output.",
    )
    benchmarkParser.add_argument(
        "-o",
        "--output",
        type=OutputDirectory("benchmark"),
        required=True,
        help="Path for folder output. Folder should not currently exist.",
    )
    benchmarkParser.add_argument(
        "-r",
        "--reference",
        type=InputFile("reference"),
        default="",
        help="Path to a reference fasta file. The reads in the input file should create a consensus that matches the reference. If no file is provided, the reference consensus sequence will be generated using all of the input reads before benchmarking.",
    )
    benchmarkParser.add_argument(
        "-c",
        "--consensusAlgorithm",
        type=ConsensusAlgorithmText(),
        default="pairwise",
        help="An option between two consensus sequence algorithms. Default is a customized algorithm that relies on pairwise alignment, which can be slow for larger sequences. Options: pairwise (default), lamassemble",
    )
    benchmarkParser.add_argument(
        "-int",
        "--intervals",
        type=BenchmarkInterval(),
        default=[10],
        help="Intervals at which benchmarking standards are set. Default is 10. For example, at default, the program will select 10 random target sequences to generate a consensus sequence, then 20 etc. Users can also provide a comma-delimited list of integers to specify specific intervals. If provided '10,25' the program will only evaluate subsamples of sizes 10 and 25.",
    )
    benchmarkParser.add_argument(
        "-iter",
        "--iterations",
        type=ConseqInt("iterations"),
        default=100,
        help="Number of iterations that occur at each interval. Default is 100. For example, at default, the program will generate a consensus sequence 100 times from randomly selecting 10 sequences, then 100 times for 20 sequences etc.",
    )
    benchmarkParser.add_argument(
        "-p",
        "--processNum",
        type=ConseqInt("processNum"),
        default=1,
        help="Number of processes to run. By default it will only use 1. If you enter a number beyond the number of processes your computer is capable of, the number of processes will automatically be set to the maximum level for your computer.",
    )
    benchmarkParser.add_argument(
        "-l",
        "--lastTrain",
        type=LastTrainFile(),
        help="Path to a last-train mat file for lamassemble. If you have already put the path to the desired file in the consensus/config.py file, this flag is unnecessary and should not be used.",
    )
    return parser


class InputDirectory:
    def __init__(self, command):
        self.allowedFileTypes = set(["fastq", "fq"])
        self.command = command

    def __call__(self, name):
        if os.path.isfile(name):
            raise argparse.ArgumentTypeError(
                "The -i or --input argument must be a directory, not a file."
            )
        if name[-1] != "/":
            name += "/"
        if not os.path.isdir(name):
            raise argparse.ArgumentTypeError(
                "The -i or --input argument must be an existing directory."
            )
        files = os.listdir(name)
        if len(files) == 0:
            raise argparse.ArgumentTypeError(
                "The -i or --input argument directory must not be empty."
            )
        for file in files:
            if file.split(".")[-1] not in self.allowedFileTypes:
                raise argparse.ArgumentTypeError(
                    f"The -i or --input argument directory must only contain fastq files (.fq or .fastq). Offending file: {file}"
                )
        if self.command == "umi":
            records = []
            for file in files:
                records.extend(list(SeqIO.parse(os.path.join(name, file), "fastq")))
            return records
        elif self.command == "cons":
            records = {}
            for file in files:
                records[os.path.join(name, file)] = list(
                    SeqIO.parse(os.path.join(name, file), "fastq")
                )
            return records


def generate_output_name(consensusAlgorithm):
    return "ConSeqUMI-" + consensusAlgorithm + time.strftime("-%Y%m%d-%H%M%S") + "/"


class OutputDirectory:
    def __init__(self, consensusAlgorithm):
        self.consensusAlgorithm = consensusAlgorithm

    def __call__(self, name):
        if os.path.isfile(name):
            raise argparse.ArgumentTypeError(
                "The -o or --output argument must be a directory, not a file."
            )
        if name[-1] != "/":
            name += "/"
        parentDirectoryPath = "/".join(name.split("/")[:-2])
        if not os.path.isdir(parentDirectoryPath):
            raise argparse.ArgumentTypeError(
                "The -o or --output argument directory requires an existing parent directory."
            )
        if not os.path.isdir(name):
            name = name[:-1] + "_" + generate_output_name(self.consensusAlgorithm)
        if os.path.isdir(name):
            name += "/" + generate_output_name(self.consensusAlgorithm)
        os.mkdir(name)
        return name


class AdapterFile:
    def __init__(self):
        self.allowedFileTypes = set(["txt"])
        self.allowedNucleotides = set([*"ACGTURYKMSWBDHVN"])

    def __call__(self, name):
        if not os.path.isfile(name):
            raise argparse.ArgumentTypeError(
                "The -a or --adapters argument must be an existing file."
            )
        if name.split(".")[-1] not in self.allowedFileTypes:
            raise argparse.ArgumentTypeError(
                "The -a or --adapters argument must be a text (.txt) file."
            )
        with open(name, "r") as adapterFile:
            adapters = [adapter.rstrip() for adapter in adapterFile.readlines()]
        if len(adapters) != 4:
            raise argparse.ArgumentTypeError(
                f"The -a or --adapters argument file must contain exactly 4 adapters. Your file contains: {len(adapters)}"
            )
        allAdapterNucleotides = []
        for adapter in adapters:
            allAdapterNucleotides.extend([*adapter])
        if len(set(allAdapterNucleotides) - self.allowedNucleotides) > 0:
            raise argparse.ArgumentTypeError(
                f"The -a or --adapters argument adapters can only contain appropriate nucleotide codes. Invalid codes: {set(allAdapterNucleotides) - self.allowedNucleotides}"
            )

        return adapters


class ConsensusAlgorithmText:
    def __init__(self):
        self.validConsensusArgorithms = set(["pairwise", "lamassemble", "medaka"])

    def __call__(self, name):
        if name not in self.validConsensusArgorithms:
            raise argparse.ArgumentTypeError(
                f"The -c or --consensusAlgorithm argument must be 'pairwise', 'lamassemble' or 'medaka'. Offending value: {name}"
            )
        if name == "medaka" and not which("medaka_consensus"):
            raise argparse.ArgumentTypeError(
                "You must install medaka in order to use the medaka consensus algorithm option."
            )
        if name == "lamassemble" and not which("lamassemble"):
            raise argparse.ArgumentTypeError(
                "You must install lamassemble in order to use the lamassemble consensus algorithm option."
            )
        return name

class LastTrainFile:

    def __call__(self, name):
        if os.path.isdir(name):
            raise argparse.ArgumentTypeError(
                f"The -l or --lastTrain argument must be a file, not a directory."
            )
        if not os.path.isfile(name):
            raise argparse.ArgumentTypeError(
                f"The -l or --lastTrain argument must be an existing file."
            )
        if name.split(".")[-1] != "mat":
            raise argparse.ArgumentTypeError(
                f"The -l or --lastTrain argument must be a .mat file."
            )
        return name


class ConseqInt:
    def __init__(self, type):
        self.minValue = 1
        if type == "minimumReads":
            self.type = "minimumReads"
            self.conciseType = "m"
        elif type == "intervals":
            self.type = "intervals"
            self.conciseType = "int"
        elif type == "iterations":
            self.type = "iterations"
            self.conciseType = "iter"
        elif type == "umiLength":
            self.minValue = 10
            self.type = "umiLength"
            self.conciseType = "u"
        elif type == "processNum":
            self.type = "processNum"
            self.conciseType = "p"

    def __call__(self, name):
        try:
            nameInt = int(name)
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"The -{self.conciseType} or --{self.type} argument must be an integer. Offending value: {name}"
            )
        if nameInt < self.minValue and not (self.type == "umiLength" and nameInt == 0):
            raise argparse.ArgumentTypeError(
                f"The -{self.conciseType} or --{self.type} argument must be greater than or equal to {self.minValue}. Offending value: {name}"
            )
        if nameInt > os.cpu_count() and self.type == "processNum":
            nameInt = None

        return nameInt


class BenchmarkInterval:
    def __call__(self, intervalInput):
        intervalInputs = intervalInput.split(",")
        conseqInt = ConseqInt("intervals")
        intervalInputs = [conseqInt(interval) for interval in intervalInputs]
        return intervalInputs


class InputFile:
    def __init__(self, type):
        if type == "input":
            self.allowedFileTypes = ["fastq", "fq"]
            self.conciseType = "i"
        elif type == "reference":
            self.allowedFileTypes = ["fasta", "fa"]
            self.conciseType = "r"
        self.type = type

    def __call__(self, name):
        if self.type == "reference" and name == "":
            return name
        if os.path.isdir(name):
            raise argparse.ArgumentTypeError(
                f"The -{self.conciseType} or --{self.type} argument must be a file, not a directory."
            )
        if not os.path.isfile(name):
            raise argparse.ArgumentTypeError(
                f"The -{self.conciseType} or --{self.type} argument must be an existing file."
            )
        if name.split(".")[-1] not in self.allowedFileTypes:
            raise argparse.ArgumentTypeError(
                f"The -{self.conciseType} or --{self.type} argument file can only be a {self.allowedFileTypes[0]} file (.{self.allowedFileTypes[1]} or .{self.allowedFileTypes[0]})."
            )
        return list(SeqIO.parse(name, self.allowedFileTypes[0]))


if __name__ == "__main__":
    main()
