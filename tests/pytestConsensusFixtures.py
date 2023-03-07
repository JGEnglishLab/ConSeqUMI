import pytest
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@pytest.fixture
def simpleInsert():
    return "C"


@pytest.fixture
def simpleString():
    return "A" * 10


@pytest.fixture
def middleInsertIndex(simpleString):
    return len(simpleString) // 2


@pytest.fixture
def consensusSequence(simpleInsert):
    random.seed(0)
    consensusSequence = "".join(random.choices("ATG", k=50))
    repeatingInsert = "ATG" * 50
    consensusSequence = (
        consensusSequence[:15] + repeatingInsert + consensusSequence[15:]
    )
    consensusSequence = (
        consensusSequence[:30] + simpleInsert * 6 + consensusSequence[36:]
    )
    return consensusSequence


def generate_read_sequences(sequence):
    insertIndex, deleteIndex, mutateIndex = 20, 30, 40
    targetSequences = []
    for i in range(1, 6):
        targetSequences.append(
            sequence[:insertIndex] + i * "C" + sequence[insertIndex:]
        )
    for i in range(1, 6):
        targetSequences.append(sequence[:deleteIndex] + sequence[deleteIndex + i :])
    for i in range(1, 5):
        targetSequences.append(
            sequence[:mutateIndex] + i * "C" + sequence[mutateIndex + i :]
        )
    return targetSequences


@pytest.fixture
def targetSequenceDifferences(simpleInsert):
    differences = []
    for i in range(1, 6):
        differences.append((20, 20, simpleInsert * i))
    for i in range(1, 6):
        differences.append((30, 30 + i, ""))
    for i in range(1, 5):
        differences.append((40, 40 + i, simpleInsert * i))
    return differences


@pytest.fixture
def targetSequences(consensusSequence):
    return generate_read_sequences(consensusSequence)


@pytest.fixture
def targetSequenceRecords(targetSequences):
    seqRecords = []
    for i in range(len(targetSequences)):
        sequence = Seq(targetSequences[i])
        id = str(i)
        phred_quality = [
            40 if targetSequences[i][j] != "C" else 10
            for j in range(len(targetSequences[i]))
        ]
        seqRecords.append(
            SeqRecord(
                sequence, id=id, letter_annotations={"phred_quality": phred_quality}
            )
        )
    return seqRecords
