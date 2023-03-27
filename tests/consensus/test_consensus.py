import pytest
import re
from Bio import SeqIO
import sys
import os

srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src/ConSeqUMI"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)

from consensus import consensus
from test_conseq import parser, consArgs, consFiles
from test_conseq import parsedConsArgs as args
from pytestConsensusFixtures import (
    consensusSequence,
    targetSequences,
    targetSequenceRecords,
    simpleInsert,
)


@pytest.fixture
def consensusAlgorithm():
    return "pairwise"


@pytest.fixture
def fileNamePattern(consensusAlgorithm):
    pattern = r"-\d{8}-\d{6}\.fasta"
    return "consensus-" + consensusAlgorithm + pattern


def test__cons__main(args, consFiles):
    consensus.main(args)
    file = os.listdir(args["output"])
    assert len(file) == 1
    consFile = args["output"] + file[0]
    expectedDescriptionStarts = [
        "Number of Target Sequences used to generate this consensus: 14",
        "Number of Target Sequences used to generate this consensus: 28",
    ]
    consensusRecords = list(SeqIO.parse(consFile, "fasta"))
    for record in consensusRecords:
        descriptionStart = len(record.id) + 1
        assert (
            record.description[
                descriptionStart : descriptionStart + len(expectedDescriptionStarts[0])
            ]
            in expectedDescriptionStarts
        )


def test__cons__main_quits_when_minimum_read_count_reached(args, consFiles):
    args["minimumReads"] = 20
    consensus.main(args)
    file = os.listdir(args["output"])
    consFile = args["output"] + file[0]
    consensusRecords = list(SeqIO.parse(consFile, "fasta"))
    assert len(consensusRecords) == 1


def test__cons__generate_file_name(consensusAlgorithm, fileNamePattern):
    fileNameOutput = consensus.generate_file_name(consensusAlgorithm)
    assert re.match(fileNamePattern, fileNameOutput)


def test__cons__generate_file_name_for_medaka_with_model():
    consensusAlgorithm = "medaka"
    pattern = r"-\d{8}-\d{6}\.fasta"
    medakaFileName = "consensus-medaka-r941_min_high_g303" + pattern
    medakaFileNameOutput = consensus.generate_file_name(consensusAlgorithm)
    assert re.match(medakaFileName, medakaFileNameOutput)


@pytest.mark.skipif(
    True,
    reason="Not sure how to test this. I would need to change the MCOMMAND config variable in-code to test separately from '_with_model'",
)
def test__cons__generate_file_name_for_medaka_with_no_model():
    consensusAlgorithm = "medaka"
    pattern = r"-\d{8}-\d{6}\.fasta"
    fileName = "consensus-medaka" + pattern
    fileNameOutput = consensus.generate_file_name(consensusAlgorithm)
    assert re.match(fileName, fileNameOutput)
