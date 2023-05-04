from ConSeqUMI.Printer import Printer
from ConSeqUMI.consensus.ConsensusContext import ConsensusContext
from ConSeqUMI.consensus.config import LCOMMAND

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import time
import argparse
import os
from concurrent.futures import ProcessPoolExecutor, Future, as_completed
import typing as T


def find_consensus_and_add_to_writing_queue(path, records, context, printer):
    printer(f" ***** {len(records)} reads: generating consensus for {path}")
    id = path.split("/")[-1]
    description = f"Number of Target Sequences used to generate this consensus: {len(records)}, File Path: {path}"
    consensusRecord = context.generate_consensus_record_from_biopython_records(records)
    consensusRecord.id = id
    consensusRecord.description = description
    return consensusRecord


def writing_to_file_from_queue(queue, consensusFilePath):
    with open(consensusFilePath, "w") as output_handle:
        while True:
            consensusRecord = queue.get()
            if consensusRecord is None:
                break
            SeqIO.write([consensusRecord], output_handle, "fasta")


def main(args):
    printer = Printer()
    context = ConsensusContext(args["consensusAlgorithm"])
    pathsSortedByLength = sorted(args["input"])
    pathsSortedByLength = sorted(
        pathsSortedByLength, key=lambda k: len(args["input"][k]), reverse=True
    )
    outputFileType = determine_output_file_type(args["consensusAlgorithm"])
    consensusFilePath = os.path.join(
        args["output"],
        context.generate_consensus_algorithm_path_header("consensus")
        + "."
        + outputFileType,
    )
    print("output folder: " + consensusFilePath)
    printer("beginning consensus sequence generation")

    consensusGenerationProcessPool: ProcessPoolExecutor = ProcessPoolExecutor(
        max_workers=args["processNum"]
    )
    futureProcesses: T.List[Future] = []

    for path in pathsSortedByLength:
        records = args["input"][path]
        if len(records) < args["minimumReads"]:
            printer(
                f"remaining files have fewer than minimum read number ({args['minimumReads']}), ending program"
            )
            break
        futureProcesses.append(
            consensusGenerationProcessPool.submit(
                find_consensus_and_add_to_writing_queue, path, records, context, printer
            )
        )

    with open(consensusFilePath, "w") as output_handle:
        for futureProcess in as_completed(futureProcesses):
            SeqIO.write([futureProcess.result()], output_handle, outputFileType)

    printer("consensus generation complete")


def determine_output_file_type(consensusAlgorithm):
    if consensusAlgorithm == "lamassemble":
        lamassembleParser = argparse.ArgumentParser(description="")
        lamassembleParser.add_argument("-f", type=str)
        lamassembleParser.add_argument("-format", type=str)
        args, unknown = lamassembleParser.parse_known_args(LCOMMAND)
        args = vars(args)
        if args["f"] or args["format"]:
            return args["f"]
    return "fasta"
