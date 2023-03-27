from ConSeqUMI.Printer import Printer
from ConSeqUMI.consensus.ConsensusContext import ConsensusContext
from ConSeqUMI.consensus.config import MCOMMAND

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import time
import argparse
import os
from multiprocessing import Process, Queue


def producer(queue, pathsSortedByLength, args, context):
    printer = Printer()
    for path in pathsSortedByLength:
        records = args["input"][path]
        if len(records) < args["minimumReads"]:
            printer(
                f"remaining files have fewer than minimum read number ({args['minimumReads']}), ending program"
            )
            break
        printer(f" ***** {len(records)} reads: generating consensus for {path}")
        id = path.split("/")[-1]
        description = f"Number of Target Sequences used to generate this consensus: {len(records)}, File Path: {path}"
        consensusSequence = (
            context.generate_consensus_sequence_from_biopython_records(records)
        )
        consensusRecord = SeqRecord(
            Seq(consensusSequence), id=id, description=description
        )
        queue.put(consensusRecord)
    queue.put(None)


def consumer(queue, consensusFilePath):
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
    consensusFilePath = os.path.join(
        args["output"],
        context.generate_consensus_algorithm_path_header("consensus") + ".fasta",
    )
    print("output folder: " + consensusFilePath)
    printer("beginning consensus sequence generation")

    queue = Queue()
    consumer_process = Process(target=consumer, args=(queue,consensusFilePath))
    consumer_process.start()
    producer_process = Process(target=producer, args=(queue,pathsSortedByLength, args, context))
    producer_process.start()
    producer_process.join()
    consumer_process.join()

    printer("consensus generation complete")




