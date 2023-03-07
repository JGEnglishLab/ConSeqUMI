import pytest
import re
import pandas as pd
import os
import sys

srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src/ConSeqUMI"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)

from consensus import benchmark
from test_conseq import parser, benchmarkArgs, benchmarkFiles
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
    pattern = r"-\d{8}-\d{6}\.csv"
    return "benchmark-" + consensusAlgorithm + pattern


def test__benchmark__main__defaults(
    parser, benchmarkArgs, benchmarkFiles, fileNamePattern
):
    benchmarkArgs += ["-iter", "1"]
    args = vars(parser.parse_args(benchmarkArgs))
    benchmark.main(args)
    outputContents = sorted(os.listdir(args["output"]))
    assert len(outputContents) == 3
    benchmarkOutput = outputContents[0]
    assert re.match(fileNamePattern, benchmarkOutput)
    assert outputContents[1] == "input.fastq"
    assert outputContents[2] == "reference.fasta"
    benchmarkOutput = args["output"] + benchmarkOutput
    benchmarkDf = pd.read_csv(benchmarkOutput)
    columns = [
        "interval",
        "iteration",
        "referenceSequence",
        "benchmarkSequence",
        "levenshteinDistance",
        "originalNumberOfSequences",
    ]
    assert list(benchmarkDf.columns) == columns
    assert len(benchmarkDf) == 2


def test__benchmark__main__works_with_provided_reference_sequence(
    parser, benchmarkArgs, benchmarkFiles, fileNamePattern
):
    benchmarkArgs += ["-r", benchmarkFiles.referenceFile.name]
    benchmarkArgs += ["-iter", "1"]
    args = vars(parser.parse_args(benchmarkArgs))
    benchmark.main(args)
    outputContents = sorted(os.listdir(args["output"]))
    assert len(outputContents) == 3
    benchmarkOutput = outputContents[0]
    assert re.match(fileNamePattern, benchmarkOutput)
    assert outputContents[1] == "input.fastq"
    assert outputContents[2] == "reference.fasta"
    benchmarkOutput = args["output"] + benchmarkOutput
    benchmarkDf = pd.read_csv(benchmarkOutput)
    columns = [
        "interval",
        "iteration",
        "referenceSequence",
        "benchmarkSequence",
        "levenshteinDistance",
        "originalNumberOfSequences",
    ]
    assert list(benchmarkDf.columns) == columns
    assert len(benchmarkDf) == 2


def test__benchmark__generate_file_name(consensusAlgorithm, fileNamePattern):
    fileNameOutput = benchmark.generate_file_name(consensusAlgorithm)
    assert re.match(fileNamePattern, fileNameOutput)
