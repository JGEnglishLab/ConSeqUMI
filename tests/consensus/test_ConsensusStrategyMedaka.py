import pytest
import sys
import os
from shutil import which
from Levenshtein import distance
import pandas as pd
import numpy as np
import re
from concurrent.futures import Future, as_completed
import typing as T

srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src/ConSeqUMI"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)
from pytestConsensusFixtures import (
    consensusSequence,
    targetSequences,
    targetSequenceRecords,
    simpleInsert,
    simpleString,
    middleInsertIndex,
    targetSequenceDifferences,
)
from consensus import ConsensusStrategyMedaka

if not which("medaka_consensus"):
    pytest.skip(
        "skipping tests that require medaka, as it has not been rendered executable on this computer.",
        allow_module_level=True,
    )


@pytest.fixture
def consensusStrategyMedaka():
    return ConsensusStrategyMedaka.ConsensusStrategyMedaka()


def test__consensus_strategy_medaka__generate_consensus_sequence_from_biopython_records(
    consensusSequence, targetSequenceRecords, consensusStrategyMedaka
):
    consensusSequenceOutput = (
        consensusStrategyMedaka.generate_consensus_record_from_biopython_records(
            targetSequenceRecords
        )
    )
    assert str(consensusSequenceOutput.seq) == consensusSequence


def test__consensus_strategy_medaka__populate_future_processes_with_benchmark_tasks(
    consensusStrategyMedaka, consensusSequence, targetSequenceRecords
):
    intervals = [10]
    iterations = 2
    numProcesses = 1
    rows = [
        ["1", "0", consensusSequence, "tempSequence", "distance", "14"],
        ["1", "1", consensusSequence, "tempSequence", "distance", "14"],
        ["10", "0", consensusSequence, "tempSequence", "distance", "14"],
        ["10", "1", consensusSequence, "tempSequence", "distance", "14"],
    ]
    futureProcesses: T.List[Future] = []
    consensusStrategyMedaka.populate_future_processes_with_benchmark_tasks(
        futureProcesses, numProcesses, consensusSequence, targetSequenceRecords, intervals, iterations
    )
    rowsOutput = []
    for futureProcess in as_completed(futureProcesses):
        row = futureProcess.result()
        rowsOutput.append(row)
    rowsOutput.sort()

    assert len(rowsOutput) == len(rows)
    for i in range(len(rows)):
        row = rows[i]
        rowOutput = rowsOutput[i]
        assert rowOutput[:3] == row[:3]
        assert distance(rowOutput[2], rowOutput[3]) == int(rowOutput[4])
        assert rowOutput[-1] == row[-1]


def test__consensus_strategy_medaka__populate_future_processes_with_benchmark_tasks__max_interval_number_is_500(
    consensusStrategyMedaka, consensusSequence, targetSequenceRecords
):
    intervals = [100]
    iterations = 1
    numberOfRecords = 605
    numProcesses = 1
    inputRecords = [targetSequenceRecords[0] for _ in range(numberOfRecords)]
    futureProcesses: T.List[Future] = []
    consensusStrategyMedaka.populate_future_processes_with_benchmark_tasks(
        futureProcesses, numProcesses, consensusSequence, inputRecords, intervals, iterations
    )
    rowsOutput = []
    for futureProcess in as_completed(futureProcesses):
        row = futureProcess.result()
        rowsOutput.append(row)
    rowsOutput.sort()
    intervalsOutput = set(pd.DataFrame(rowsOutput).iloc[:, 0])
    intervals = set(["1", "100", "200", "300", "400", "500"])
    assert intervalsOutput == intervals


def test__consensus_strategy_medaka__populate_future_processes_with_benchmark_tasks__customized_intervals_also_works(
    consensusStrategyMedaka, consensusSequence, targetSequenceRecords
):
    intervals = [7, 10, 12]
    iterations = 2
    numProcesses = 1
    rows = [
        ["10", "0", consensusSequence, "tempSequence", "distance", "14"],
        ["10", "1", consensusSequence, "tempSequence", "distance", "14"],
        ["12", "0", consensusSequence, "tempSequence", "distance", "14"],
        ["12", "1", consensusSequence, "tempSequence", "distance", "14"],
        ["7", "0", consensusSequence, "tempSequence", "distance", "14"],
        ["7", "1", consensusSequence, "tempSequence", "distance", "14"],
    ]
    futureProcesses: T.List[Future] = []
    consensusStrategyMedaka.populate_future_processes_with_benchmark_tasks(
        futureProcesses, numProcesses, consensusSequence, targetSequenceRecords, intervals, iterations
    )
    rowsOutput = []
    for futureProcess in as_completed(futureProcesses):
        row = futureProcess.result()
        rowsOutput.append(row)
    rowsOutput.sort()

    assert len(rowsOutput) == len(rows)
    for i in range(len(rows)):
        row = rows[i]
        rowOutput = rowsOutput[i]
        assert rowOutput[:3] == row[:3]
        assert distance(rowOutput[2], rowOutput[3]) == int(rowOutput[4])
        assert rowOutput[-1] == row[-1]


def test__consensus_strategy_medaka__generate_consensus_algorithm_path_header__for_medaka_with_model(
    consensusStrategyMedaka,
):
    processName = "consensus"
    medakaFileName = processName + r"-medaka-r941_min_high_g303-\d{8}-\d{6}"
    medakaFileNameOutput = (
        consensusStrategyMedaka.generate_consensus_algorithm_path_header(processName)
    )
    assert re.match(medakaFileName, medakaFileNameOutput)


@pytest.mark.skipif(
    True,
    reason="Not sure how to test this. I would need to change the MCOMMAND config variable in-code to test separately from '_with_model'",
)
def test__consensus_strategy_medaka__generate_consensus_algorithm_path_header__for_medaka_with_no_model(
    consensusStrategyMedaka,
):
    processName = "consensus"
    medakaFileName = processName + r"-medaka-\d{8}-\d{6}"
    medakaFileNameOutput = (
        consensusStrategyMedaka.generate_consensus_algorithm_path_header(processName)
    )
    assert re.match(medakaFileName, medakaFileNameOutput)
