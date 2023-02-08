import pytest
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)
from pytestConsensusFixtures import consensusSequence, targetSequences, targetSequenceRecords, simpleInsert, simpleString, middleInsertIndex, targetSequenceDifferences
from consensus import ConsensusStrategyMedaka

@pytest.fixture
def consensusStrategyMedaka():
    return ConsensusStrategyMedaka.ConsensusStrategyMedaka()
