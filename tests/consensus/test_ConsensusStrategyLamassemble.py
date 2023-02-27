import pytest
from Levenshtein import distance
import sys
import os
from shutil import which

srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)
from pytestConsensusFixtures import consensusSequence, targetSequences, targetSequenceRecords, simpleInsert, simpleString, middleInsertIndex, targetSequenceDifferences
from consensus import ConsensusStrategyLamassemble

if not which("lamassemble"):
    pytest.skip("skipping tests that require lamassemble, as it has not been rendered executable on this computer.", allow_module_level=True)

@pytest.fixture
def consensusStrategyLamassemble():
    return ConsensusStrategyLamassemble.ConsensusStrategyLamassemble()

def test__consensus_strategy_lamassemble__generate_consensus_sequence_from_biopython_records(consensusSequence, targetSequenceRecords, consensusStrategyLamassemble):
    consensusSequenceOutput = consensusStrategyLamassemble.generate_consensus_sequence_from_biopython_records(targetSequenceRecords)
    assert consensusSequenceOutput == consensusSequence

def test__consensus_strategy_lamassemble__benchmark_sequence_generator(consensusStrategyLamassemble, consensusSequence, targetSequenceRecords):
    intervals = 10
    iterations = 2
    rows = [
        ["1","0",consensusSequence,"tempSequence","distance","14"],
        ["1","1",consensusSequence,"tempSequence","distance","14"],
        ["10","0",consensusSequence,"tempSequence","distance","14"],
        ["10","1",consensusSequence,"tempSequence","distance","14"],
    ]
    rowsOutput = [row for row in consensusStrategyLamassemble.benchmark_sequence_generator(consensusSequence, targetSequenceRecords, intervals, iterations)]

    assert len(rowsOutput) == len(rows)
    for i in range(len(rows)):
        row = rows[i]
        rowOutput = rowsOutput[i]
        assert rowOutput[:3] == row[:3]
        assert distance(rowOutput[2],rowOutput[3]) == int(rowOutput[4])
        assert rowOutput[-1] == row[-1]

def test__consensus_strategy_lamassemble__benchmark_sequence_generator__max_number_of_intervals_is_100(consensusStrategyLamassemble, consensusSequence, targetSequenceRecords):
    intervals = 1
    iterations = 1
    numberOfRecords = 101
    inputRecords = [targetSequenceRecords[0] for _ in range(numberOfRecords)]
    rowsOutput = [row for row in consensusStrategyLamassemble.benchmark_sequence_generator(consensusSequence, inputRecords, intervals, iterations)]
    assert len(rowsOutput) == 100
