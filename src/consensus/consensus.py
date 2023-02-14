from Printer import Printer
from consensus.ConsensusContext import ConsensusContext
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import time

def main(args):
    printer = Printer()
    context = ConsensusContext(args["consensusAlgorithm"])
    pathsSortedByLength = sorted(args["input"])
    pathsSortedByLength = sorted(pathsSortedByLength, key=lambda k: len(args["input"][k]), reverse=True)
    consensusFilePath = "/".join(pathsSortedByLength[0].split("/")[:-2]) + "/" + generate_file_name(args["consensusAlgorithm"])
    print(consensusFilePath)
    printer("beginning consensus sequence generation")
    with open(consensusFilePath, "w") as output_handle:
        for path in pathsSortedByLength:
            records = args["input"][path]
            if len(records) < args["minimumReads"]:
                printer(f"remaining files have fewer than minimum read number ({args['minimumReads']}), ending program")
                break
            printer(f" ***** {len(records)} reads: generating consensus for {path}")
            id = path.split("/")[-1]
            description = f"Number of Target Sequences used to generate this consensus: {len(records)}, File Path: {path}"
            consensusSequence = context.generate_consensus_sequence_from_biopython_records(records)
            consensusRecord = SeqRecord(Seq(consensusSequence), id=id, description=description)
            SeqIO.write([consensusRecord], output_handle, "fasta")
    printer("consensus generation complete")


def generate_file_name(consensusAlgorithm):
    return "consensus-" + consensusAlgorithm + time.strftime("-%Y%m%d-%H%M%S.fasta")
