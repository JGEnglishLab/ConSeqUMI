import pytest
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@pytest.fixture
def consensusSequence():
    random.seed(0)
    consensusSequence = "".join(random.choices("ATGC", k=30))
    repeatingInsert = "ATGC" * 5
    consensusSequence = consensusSequence[:15] + repeatingInsert + consensusSequence[15:]
    return consensusSequence

def generate_read_sequences(sequence):
    readSequences = []
    insertIndex = 20
    for i in range(1, 6):
        readSequences.append(sequence[:insertIndex] + i*"A" + sequence[insertIndex:])
    deleteIndex = 30
    for i in range(1, 6):
        readSequences.append(sequence[:deleteIndex] + sequence[deleteIndex + i:])
    mutationIndex = 40
    for i in range(5):
        readSequences.append(sequence[:mutationIndex + i] + "R" + sequence[mutationIndex + i + 1:])
    return readSequences

@pytest.fixture
def readSequences(consensusSequence):
    return generate_read_sequences(consensusSequence)

@pytest.fixture
def readSequenceRecords(readSequences):
    return [SeqRecord(Seq(readSequences[i]), id=str(i)) for i in range(len(readSequences))]
