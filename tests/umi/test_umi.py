import pytest
import random
from umi import umi
from pytestConsensusFixtures import consensusSequence, targetSequences, targetSequenceRecords, simpleInsert, simpleString, middleInsertIndex, targetSequenceDifferences




@pytest.fixture
def readRecords():
    consensusSequence = "A"
