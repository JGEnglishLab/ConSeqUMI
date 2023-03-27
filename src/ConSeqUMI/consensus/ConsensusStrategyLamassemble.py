from ConSeqUMI.consensus.ConsensusStrategy import ConsensusStrategy
import subprocess
from Bio import SeqIO
from io import StringIO
from tempfile import NamedTemporaryFile
from ConSeqUMI.consensus.config import LCOMMAND


class ConsensusStrategyLamassemble(ConsensusStrategy):
    def generate_consensus_algorithm_path_header_insert(self) -> str:
        return "lamassemble"

    def generate_consensus_sequence_from_biopython_records(
        self, binRecords: list
    ) -> str:
        inputFile = NamedTemporaryFile(
            prefix="conseq_lamassemble_delete_", suffix=".fastq"
        )
        with open(inputFile.name, "w") as output_handle:
            SeqIO.write(binRecords, output_handle, "fastq")

        processCommands = LCOMMAND[:] + [inputFile.name]
        child = subprocess.Popen(
            processCommands,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        # child.stdin.write(stdin.encode())
        child_out = child.communicate()[0].decode("utf8")
        seq_ali = list(SeqIO.parse(StringIO(child_out), "fasta"))
        child.stdin.close()

        return str(seq_ali[0].seq).upper()
