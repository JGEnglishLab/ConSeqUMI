from ConSeqUMI.Printer import Printer
from ConSeqUMI.consensus.ConsensusContext import ConsensusContext
import time
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
from concurrent.futures import Future, as_completed
import typing as T

def writing_to_file_from_queue(queue, benchmarkOutputFile):
    with open(benchmarkOutputFile, "w") as file:
        while True:
            row = queue.get()
            if row is None:
                break
            file.write(",".join(row) + os.linesep)


def main(args):
    context = ConsensusContext(args["consensusAlgorithm"])
    inputFile = os.path.join(args["output"], "input.fastq")
    benchmarkOutputFile = os.path.join(
        args["output"],
        context.generate_consensus_algorithm_path_header("benchmark") + ".csv",
    )
    print("output folder: " + benchmarkOutputFile)
    referenceFile = os.path.join(args["output"], "reference.fasta")
    printer = Printer()
    print(args["output"])
    printer(f"total number of input reads: {len(args['input'])}")
    columns = [
        "interval",
        "iteration",
        "referenceSequence",
        "benchmarkSequence",
        "levenshteinDistance",
        "originalNumberOfSequences",
    ]
    if args["reference"]:
        referenceSequence = str(args["reference"][0].seq)
    else:
        printer("no reference sequence provided. Generating reference sequence")
        referenceSequence = context.generate_consensus_sequence_from_biopython_records(
            args["input"]
        )

    printer("writing input and reference sequence values to file for future reference")
    with open(inputFile, "w") as output_handle:
        SeqIO.write(args["input"], output_handle, "fastq")
    with open(referenceFile, "w") as output_handle:
        SeqIO.write(
            SeqRecord(Seq(referenceSequence), id="reference"), output_handle, "fasta"
        )

    printer("beginning benchmark process")

    futureProcesses: T.List[Future] = []
    context.populate_future_processes_with_benchmark_tasks(futureProcesses, referenceSequence, args["input"], args["intervals"], args["iterations"])
    priorIntervals = set()
    with open(benchmarkOutputFile, "w") as file:
        file.write(",".join(columns) + os.linesep)
        for futureProcess in as_completed(futureProcesses):
            row = futureProcess.result()
            if row[0] not in priorIntervals:
                priorIntervals.add(row[0])
                printer(
                    f"benchmarking interval: {row[0]} ({args['iterations']} iterations)"
                )
            file.write(",".join(row) + os.linesep)