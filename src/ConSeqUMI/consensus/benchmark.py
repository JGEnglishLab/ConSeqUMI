from ConSeqUMI.Printer import Printer
from ConSeqUMI.consensus.ConsensusContext import ConsensusContext
import time
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def main(args):
    inputFile = args["output"] + "input.fastq"
    benchmarkOutputFile = args["output"] + generate_file_name(args["consensusAlgorithm"])
    print("output folder: " + benchmarkOutputFile)
    referenceFile = args["output"] + "reference.fasta"
    printer = Printer()
    context = ConsensusContext(args["consensusAlgorithm"])
    print(args["output"])
    printer(f"total number of input reads: {len(args['input'])}")
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
        printer("no reference sequence provided. Generating reference sequence")
        referenceSequence = context.generate_consensus_sequence_from_biopython_records(args["input"])

    printer("writing input and reference sequence values to file for future reference")
    with open(inputFile, "w") as output_handle:
        SeqIO.write(args["input"], output_handle, "fastq")
    with open(referenceFile, "w") as output_handle:
        SeqIO.write(SeqRecord(Seq(referenceSequence), id="reference"), output_handle, "fasta")

    printer("beginning benchmark process")
    with open(benchmarkOutputFile, "w") as file:
        file.write(",".join(columns) + os.linesep)
        previousInterval = 0
        for outputValue in context.benchmark_sequence_generator(referenceSequence, args["input"], args["intervals"], args["iterations"]):
            currentInterval = int(outputValue[0])
            if currentInterval > previousInterval:
                printer(f"benchmarking interval: {currentInterval} ({args['iterations']} iterations)")
                previousInterval = currentInterval
            file.write(",".join(outputValue) + os.linesep)

def generate_file_name(consensusAlgorithm):
    return "benchmark-" + consensusAlgorithm + time.strftime("-%Y%m%d-%H%M%S.csv")