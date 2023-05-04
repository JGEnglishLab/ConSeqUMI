from ConSeqUMI.consensus.ConsensusStrategy import ConsensusStrategy
from ConSeqUMI.consensus.ReferenceConsensusGenerator import ReferenceConsensusGenerator
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from io import StringIO
from tempfile import TemporaryDirectory, NamedTemporaryFile
from os.path import exists
from ConSeqUMI.consensus.config import MCOMMAND
import argparse


class ConsensusStrategyMedaka(ConsensusStrategy):
    def generate_consensus_algorithm_path_header_insert(self) -> str:
        consensusAlgorithmInsert = "medaka"
        medakaParser = argparse.ArgumentParser(description="")
        medakaParser.add_argument("-m", type=str)
        args, unknown = medakaParser.parse_known_args(MCOMMAND)
        args = vars(args)
        if args["m"]:
            consensusAlgorithmInsert += "-" + args["m"]
        return consensusAlgorithmInsert

    def generate_consensus_record_from_biopython_records(self, binRecords: list) -> str:
        inputFile = NamedTemporaryFile(prefix="conseq_medaka_delete_", suffix=".fastq")
        outputDir = TemporaryDirectory(prefix="conseq_medaka_delete_")
        draftFile = NamedTemporaryFile(prefix="conseq_medaka_delete_", suffix=".fasta")
        with open(inputFile.name, "w") as output_handle:
            SeqIO.write(binRecords, output_handle, "fastq")
        binSequences = [str(record.seq) for record in binRecords]
        referenceConsensusGenerator = ReferenceConsensusGenerator()
        referenceSequence = referenceConsensusGenerator.generate_consensus_sequence(
            binSequences
        )
        consensusRecords = [SeqRecord(Seq(referenceSequence), id="medaka_draft")]
        while len(consensusRecords) < 5:
            with open(draftFile.name, "w") as output_handle:
                SeqIO.write(
                    [consensusRecords[-1]],
                    output_handle,
                    "fasta",
                )
            processCommands = MCOMMAND[:]
            processCommands += [
                "-i",
                inputFile.name,
                "-d",
                draftFile.name,
                "-o",
                outputDir.name,
            ]
            child = subprocess.Popen(
                processCommands,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
            )
            child.stdin.close()
            if exists(outputDir.name + "/consensus.fasta"):
                medakaOutputRecords = [
                    record
                    for record in SeqIO.parse(
                        outputDir.name + "/consensus.fasta", "fasta"
                    )
                ]
                if len(medakaOutputRecords) == 0:
                    return consensusRecords[-1]
                consensusRecords.append(medakaOutputRecords[0])
                if str(consensusRecords[-1].seq) == str(consensusRecords[-2].seq):
                    return consensusRecords[-1]
            else:
                return consensusRecords[-1]
        return consensusRecords[-1]
