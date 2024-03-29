import pytest
import random
from Levenshtein import distance
import sys
import os
import numpy as np
import pandas as pd
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
from consensus import ConsensusStrategyPairwise


@pytest.fixture
def consensusStrategyPairwise():
    return ConsensusStrategyPairwise.ConsensusStrategyPairwise()


def test__consensus_strategy_pairwise__initialization(consensusStrategyPairwise):
    assert consensusStrategyPairwise.aligner.mismatch_score == -1
    assert consensusStrategyPairwise.aligner.open_gap_score == -1
    assert consensusStrategyPairwise.aligner.extend_gap_score == -0.5


@pytest.fixture
def stretchLength():
    return 4


@pytest.fixture
def simpleStringWithThreeInserts(
    simpleInsert, simpleString, middleInsertIndex, stretchLength
):
    simpleStringWithThreeInserts = (
        simpleInsert * stretchLength
        + simpleString[:middleInsertIndex]
        + simpleInsert
        + simpleString[middleInsertIndex:]
        + simpleInsert * stretchLength
    )
    return simpleStringWithThreeInserts


def test__consensus_strategy_pairwise__find_pairwise_score_and_all_differences_between_two_sequences__finds_all_inserts(
    simpleInsert,
    simpleString,
    middleInsertIndex,
    stretchLength,
    simpleStringWithThreeInserts,
    consensusStrategyPairwise,
):
    insertionDifferences = [
        (0, 0, simpleInsert * stretchLength),
        (middleInsertIndex, middleInsertIndex, simpleInsert),
        (len(simpleString), len(simpleString), simpleInsert * stretchLength),
    ]
    (
        _,
        insertionDifferencesOutput,
    ) = consensusStrategyPairwise.find_pairwise_score_and_all_differences_between_two_sequences(
        simpleString, simpleStringWithThreeInserts
    )
    assert insertionDifferencesOutput == insertionDifferences


def test__consensus_strategy_pairwise__find_pairwise_score_and_all_differences_between_two_sequences__finds_all_deletions(
    simpleInsert,
    simpleString,
    middleInsertIndex,
    stretchLength,
    simpleStringWithThreeInserts,
    consensusStrategyPairwise,
):
    deletionDifferences = [
        (0, stretchLength, ""),
        (
            middleInsertIndex + stretchLength,
            middleInsertIndex + stretchLength + len(simpleInsert),
            "",
        ),
        (
            len(simpleStringWithThreeInserts) - stretchLength,
            len(simpleStringWithThreeInserts),
            "",
        ),
    ]
    (
        _,
        deletionDifferencesOutput,
    ) = consensusStrategyPairwise.find_pairwise_score_and_all_differences_between_two_sequences(
        simpleStringWithThreeInserts, simpleString
    )
    assert deletionDifferencesOutput == deletionDifferences


def test__consensus_strategy_pairwise__find_pairwise_score_and_all_differences_between_two_sequences__finds_all_mutations(
    simpleInsert,
    simpleString,
    middleInsertIndex,
    stretchLength,
    simpleStringWithThreeInserts,
    consensusStrategyPairwise,
):
    simpleStringWithThreeMutations = (
        simpleInsert
        + simpleString[1:middleInsertIndex]
        + simpleInsert
        + simpleString[middleInsertIndex + 1 : -1]
        + simpleInsert
    )
    mutationDifferences = [
        (0, 1, simpleInsert),
        (middleInsertIndex, middleInsertIndex + len(simpleInsert), simpleInsert),
        (len(simpleString) - 1, len(simpleString), simpleInsert),
    ]
    (
        _,
        mutationDifferencesOutput,
    ) = consensusStrategyPairwise.find_pairwise_score_and_all_differences_between_two_sequences(
        simpleString, simpleStringWithThreeMutations
    )
    assert mutationDifferencesOutput == mutationDifferences


def test__consensus_strategy_pairwise__find_pairwise_score_and_all_differences_between_two_sequences__finds_all_deletions_after_insertion(
    simpleInsert,
    simpleString,
    middleInsertIndex,
    stretchLength,
    simpleStringWithThreeInserts,
    consensusStrategyPairwise,
):
    simpleStringWithLongerFrontInsert = (
        simpleInsert * (stretchLength + 1) + simpleString
    )
    deletionDifferences = [
        (0, 0, simpleInsert),
        (
            middleInsertIndex + stretchLength,
            middleInsertIndex + stretchLength + len(simpleInsert),
            "",
        ),
        (
            len(simpleStringWithThreeInserts) - stretchLength,
            len(simpleStringWithThreeInserts),
            "",
        ),
    ]
    (
        _,
        deletionDifferencesOutput,
    ) = consensusStrategyPairwise.find_pairwise_score_and_all_differences_between_two_sequences(
        simpleStringWithThreeInserts, simpleStringWithLongerFrontInsert
    )
    assert deletionDifferencesOutput == deletionDifferences


def calculate_average_pairwise_alignment_score_for_tests(consensusSequence):
    baseScore = len(consensusSequence)
    initialErrorScore = baseScore - 1
    allScores = []
    insertionScores = [baseScore - i * 0.5 - 1 for i in range(5)]
    deletionScores = [baseScore - i * 1.5 - 2 for i in range(5)]
    mutationScores = [baseScore - i * 2 for i in range(1, 5)]
    allScores.extend(insertionScores)
    allScores.extend(deletionScores)
    allScores.extend(mutationScores)
    return sum(allScores) / len(allScores)


def test__consensus_strategy_pairwise__find_average_pairwise_alignment_score_and_all_differences_between_candidate_sequence_and_binned_sequences__finds_average_score(
    consensusSequence, targetSequences, consensusStrategyPairwise
):
    averagePairwiseAlignmentScore = (
        calculate_average_pairwise_alignment_score_for_tests(consensusSequence)
    )
    (
        averagePairwiseAlignmentScoreOutput,
        _,
    ) = consensusStrategyPairwise.find_average_pairwise_alignment_score_and_all_differences_between_candidate_sequence_and_binned_sequences(
        consensusSequence, targetSequences
    )
    assert averagePairwiseAlignmentScore == averagePairwiseAlignmentScoreOutput


def test__consensus_strategy_pairwise__find_average_pairwise_alignment_score_and_all_differences_between_candidate_sequence_and_binned_sequences__finds_all_differences(
    consensusSequence,
    targetSequences,
    consensusStrategyPairwise,
    targetSequenceDifferences,
):
    (
        _,
        targetSequenceDifferencesOutput,
    ) = consensusStrategyPairwise.find_average_pairwise_alignment_score_and_all_differences_between_candidate_sequence_and_binned_sequences(
        consensusSequence, targetSequences
    )
    assert targetSequenceDifferencesOutput == targetSequenceDifferences


def test__consensus_strategy_pairwise__generate_consensus_record_from_biopython_records(
    consensusSequence, targetSequences, consensusStrategyPairwise, targetSequenceRecords
):
    halfwayIndex = len(consensusSequence) // 2
    referenceSequence = (
        consensusSequence[:halfwayIndex] + consensusSequence[halfwayIndex + 4 :]
    )
    consensusSequenceOutput = (
        consensusStrategyPairwise.generate_consensus_record_from_biopython_records(
            targetSequenceRecords
        )
    )
    assert str(consensusSequenceOutput.seq) != referenceSequence
    assert str(consensusSequenceOutput.seq) == consensusSequence


def test__consensus_strategy_pairwise__generate_consensus_sequence_from_biopython_records__works_when_all_target_sequences_are_the_same(
    consensusStrategyPairwise, targetSequenceRecords
):
    identicalTargetSequenceRecords = [targetSequenceRecords[0] for _ in range(10)]
    expectedSequence = str(targetSequenceRecords[0].seq)
    expectedSequenceOutput = (
        consensusStrategyPairwise.generate_consensus_record_from_biopython_records(
            identicalTargetSequenceRecords
        )
    )
    assert str(expectedSequenceOutput.seq) == expectedSequence


def test__consensus_strategy_pairwise__populate_future_processes_with_benchmark_tasks(
    consensusStrategyPairwise, consensusSequence, targetSequenceRecords
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
    consensusStrategyPairwise.populate_future_processes_with_benchmark_tasks(
        futureProcesses,
        numProcesses,
        consensusSequence,
        targetSequenceRecords,
        intervals,
        iterations,
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


def test__consensus_strategy_pairwise__populate_future_processes_with_benchmark_tasks__max_interval_number_is_500(
    consensusStrategyPairwise, consensusSequence, targetSequenceRecords
):
    intervals = [100]
    iterations = 1
    numberOfRecords = 605
    numProcesses = 1
    inputRecords = [targetSequenceRecords[0] for _ in range(numberOfRecords)]
    futureProcesses: T.List[Future] = []
    consensusStrategyPairwise.populate_future_processes_with_benchmark_tasks(
        futureProcesses,
        numProcesses,
        consensusSequence,
        inputRecords,
        intervals,
        iterations,
    )
    rowsOutput = []
    for futureProcess in as_completed(futureProcesses):
        row = futureProcess.result()
        rowsOutput.append(row)
    rowsOutput.sort()
    intervalsOutput = set(pd.DataFrame(rowsOutput).iloc[:, 0])
    intervals = set(["1", "100", "200", "300", "400", "500"])
    assert intervalsOutput == intervals


def test__consensus_strategy_pairwise__populate_future_processes_with_benchmark_tasks__customized_intervals_also_works(
    consensusStrategyPairwise, consensusSequence, targetSequenceRecords
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
    consensusStrategyPairwise.populate_future_processes_with_benchmark_tasks(
        futureProcesses,
        numProcesses,
        consensusSequence,
        targetSequenceRecords,
        intervals,
        iterations,
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


def test__consensus_strategy_pairwise__generate_consensus_algorithm_path_header(
    consensusStrategyPairwise,
):
    processName = "consensus"
    medakaFileName = processName + r"-pairwise-\d{8}-\d{6}"
    medakaFileNameOutput = (
        consensusStrategyPairwise.generate_consensus_algorithm_path_header(processName)
    )
    assert re.match(medakaFileName, medakaFileNameOutput)
