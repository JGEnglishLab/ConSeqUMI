import pytest
import random
from umi import umi

def make_string_to_fastq_record(string, id):
    sequence = Seq(string)
    id = str(id)
    phred_quality = [40 for j in range(len(string))]
    return SeqRecord(sequence, id=id, letter_annotations={"phred_quality":phred_quality})
