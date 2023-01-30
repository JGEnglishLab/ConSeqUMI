import pytest
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

@pytest.fixture
def simpleInsert():
    return "R"

@pytest.fixture
def simpleString():
    return "A"*10

@pytest.fixture
def middleInsertIndex(simpleString):
    return len(simpleString) // 2

@pytest.fixture
def consensusSequence(simpleInsert):
    random.seed(0)
    consensusSequence = "".join(random.choices("ATGC", k=30))
    repeatingInsert = "ATGC" * 10
    consensusSequence = consensusSequence[:15] + repeatingInsert + consensusSequence[15:]
    consensusSequence = consensusSequence[:30] + simpleInsert*6 + consensusSequence[36:]
    return consensusSequence

def generate_read_sequences(sequence):
    insertIndex, deleteIndex, mutateIndex = 20, 30, 40
    readSequences = []
    for i in range(1, 6):
        readSequences.append(sequence[:insertIndex] + i*"R" + sequence[insertIndex:])
    for i in range(1, 6):
        readSequences.append(sequence[:deleteIndex] + sequence[deleteIndex + i:])
    for i in range(1, 5):
        readSequences.append(sequence[:mutateIndex] + i*"R" + sequence[mutateIndex + i:])
    return readSequences

@pytest.fixture
def readSequenceDifferences(simpleInsert):
    differences = []
    for i in range(1, 6):
        differences.append((20, 20, simpleInsert*i))
    for i in range(1, 6):
        differences.append((30, 30+i, ""))
    for i in range(1, 5):
        differences.append((40, 40+i, simpleInsert*i))
    return differences


@pytest.fixture
def readSequences(consensusSequence):
    return generate_read_sequences(consensusSequence)

@pytest.fixture
def readSequenceRecords(readSequences):
    return [SeqRecord(Seq(readSequences[i]), id=str(i)) for i in range(len(readSequences))]
