import sys
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
from consensus.ConsensusStrategy import ConsensusStrategy, convert_biopython_fastq_record_to_string

def test__consensus_strategy__convert_biopython_fastq_record_to_string():
    sequence = "AAAAA"
    letter_annotations = {"phred_quality":[16,16,16,16,16]}
    letterAnnotationsDecoded = ["1","1","1","1","1"]
    id = "test_id"
    fastqRecordString = f"@{id} <unknown description>\n{sequence}\n+\n{''.join(letterAnnotationsDecoded)}\n"
    fastqRecord = SeqRecord(Seq(sequence), id=id, letter_annotations=letter_annotations)
    fastqRecordStringOutput = convert_biopython_fastq_record_to_string(fastqRecord)
    assert fastqRecordStringOutput == fastqRecordString
