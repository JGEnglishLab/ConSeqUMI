from ConSeqUMI.consensus.ConsensusStrategy import ConsensusStrategy
import subprocess
from Bio import SeqIO
from io import StringIO
import tempfile

class ConsensusStrategyMedaka(ConsensusStrategy):
    def generate_consensus_sequence_from_biopython_records(self, binRecords: list) -> str:
        pass
