from consensus.ConsensusStrategy import ConsensusStrategy
import subprocess
from Bio import SeqIO
from io import StringIO
import tempfile

class ConsensusStrategyLamassemble(ConsensusStrategy):
    def generate_consensus_sequence_from_biopython_records(self, binRecords: list) -> str:

        lamassembleCommandLine = "lamassemble dependencies_download/promethion.mat --end -g60 -m 40"
        LCOMMAND = lamassembleCommandLine.split()
        with tempfile.NamedTemporaryFile(suffix=".fastq") as fp:
            with open(fp.name, "w") as output_handle:
                SeqIO.write(binRecords, output_handle, "fastq")
            processCommands = LCOMMAND[:] + [fp.name]
            child = subprocess.Popen(
                processCommands,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
            )
            #child.stdin.write(stdin.encode())
            child_out = child.communicate()[0].decode("utf8")
        seq_ali = list(SeqIO.parse(StringIO(child_out), "fasta"))
        child.stdin.close()

        return str(seq_ali[0].seq).upper()
