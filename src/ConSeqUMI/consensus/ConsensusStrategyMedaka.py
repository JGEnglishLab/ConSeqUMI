from ConSeqUMI.consensus.ConsensusStrategy import ConsensusStrategy
from ConSeqUMI.consensus.ReferenceConsensusGenerator import ReferenceConsensusGenerator
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from io import StringIO
from tempfile import TemporaryDirectory, NamedTemporaryFile
from os.path import exists

# from ConSeqUMI.consensus.config import MCOMMAND


class ConsensusStrategyMedaka(ConsensusStrategy):
    def generate_consensus_sequence_from_biopython_records(
        self, binRecords: list
    ) -> str:
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
        consensusSequences = [referenceSequence]
        while len(consensusSequences) < 5:
            with open(draftFile.name, "w") as output_handle:
                SeqIO.write(
                    [SeqRecord(Seq(consensusSequences[-1]), id="medaka_draft")],
                    output_handle,
                    "fasta",
                )
            processCommands = [
                "medaka_consensus",
                "-i",
                inputFile.name,
                "-d",
                draftFile.name,
                "-o",
                outputDir.name,
                "-f",
            ]
            child = subprocess.Popen(
                processCommands,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
            )
            child.stdin.close()
            if exists(outputDir.name + "/consensus.fasta"):
                consensusRecords = [
                    record
                    for record in SeqIO.parse(
                        outputDir.name + "/consensus.fasta", "fasta"
                    )
                ]
                if len(consensusRecords) == 0:
                    return consensusSequences[-1]
                consensusSequences.append(str(consensusRecords[0].seq))
                if consensusSequences[-1] == consensusSequences[-2]:
                    return consensusSequences[-1]
            else:
                return consensusSequences[-1]
        return consensusSequences[-1]
