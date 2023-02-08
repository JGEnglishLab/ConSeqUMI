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
from consensus import ConsensusStrategyLamassemble

@pytest.fixture
def consensusStrategyLamassemble():
    return ConsensusStrategyLamassemble.ConsensusStrategyLamassemble()

def test__consensus_strategy_lamassemble__generate_consensus_sequence_from_biopython_records(consensusSequence, targetSequenceRecords, consensusStrategyLamassemble):
    consensusSequenceOutput = consensusStrategyLamassemble.generate_consensus_sequence_from_biopython_records(targetSequenceRecords)
    assert consensusSequenceOutput == consensusSequence
