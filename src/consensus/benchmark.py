from Printer import Printer
from consensus.ConsensusContext import ConsensusContext
import time
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def main(args):
    inputFile = args["output"] + "input.fastq"
    benchmarkOutputFile = args["output"] + generate_file_name(args["consensusAlgorithm"])
    referenceFile = args["output"] + "reference.fasta"
    printer = Printer()
    context = ConsensusContext(args["consensusAlgorithm"])
    columns = [
        "interval",
        "iteration",
        "referenceSequence",
        "benchmarkSequence",
        "levenshteinDistance",
        "originalNumberOfSequences"
    ]
    if args["reference"]:
        referenceSequence = str(args["reference"][0].seq)
    else:
        referenceSequence = context.generate_consensus_sequence_from_biopython_records(args["input"])

    with open(inputFile, "w") as output_handle:
        SeqIO.write(args["input"], output_handle, "fastq")
    with open(referenceFile, "w") as output_handle:
        SeqIO.write(SeqRecord(Seq(referenceSequence), id="reference"), output_handle, "fasta")

    with open(benchmarkOutputFile, "w") as file:
        file.write(",".join(columns) + os.linesep)
        for outputValue in context.benchmark_sequence_generator(referenceSequence, args["input"], args["intervals"], args["iterations"]):
            file.write(",".join(outputValue) + os.linesep)

def generate_file_name(consensusAlgorithm):
    return "benchmark-" + consensusAlgorithm + time.strftime("-%Y%m%d-%H%M%S.csv")