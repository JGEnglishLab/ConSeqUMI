import pytest
import re
from Bio import SeqIO
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)

from consensus import consensus
from test_conseq import parser, consArgs, consFiles
from test_conseq import parsedConsArgs as args
from pytestConsensusFixtures import consensusSequence, targetSequences, targetSequenceRecords, simpleInsert

@pytest.fixture
def consensusAlgorithm():
    return "pairwise"

@pytest.fixture
def fileNamePattern(consensusAlgorithm):
    pattern = r"-\d{8}-\d{6}\.fasta"
    return "consensus-" + consensusAlgorithm + pattern

def test__cons__main(args, consFiles):
    consensus.main(args)
    fileAndDir = os.listdir(consFiles.parentDir.name)
    assert len(fileAndDir) == 2
    if fileAndDir[0].split(".")[-1] == "fasta": consFile = consFiles.parentDir.name + "/" + fileAndDir[0]
    else: consFile = consFiles.parentDir.name + "/" + fileAndDir[1]
    expectedDescriptionStarts = [
        "Number of Target Sequences used to generate this consensus: 14",
        "Number of Target Sequences used to generate this consensus: 28",
    ]
    consensusRecords = list(SeqIO.parse(consFile, "fasta"))
    for record in consensusRecords:
        descriptionStart = len(record.id)+1
        assert record.description[descriptionStart:descriptionStart + len(expectedDescriptionStarts[0])] in expectedDescriptionStarts

def test__cons__generate_file_name(consensusAlgorithm, fileNamePattern):
    fileNameOutput = consensus.generate_file_name(consensusAlgorithm)
    assert re.match(fileNamePattern, fileNameOutput)



