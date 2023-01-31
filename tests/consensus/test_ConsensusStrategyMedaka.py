from pytestConsensusFixtures import consensusSequence, readSequences, readSequenceRecords, simpleInsert, simpleString, middleInsertIndex, readSequenceDifferences
import pytest
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
from consensus import ConsensusStrategyMedaka

@pytest.fixture
def consensusStrategyMedaka():
    return ConsensusStrategyMedaka.ConsensusStrategyMedaka()
