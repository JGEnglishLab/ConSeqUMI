import pytest
import random
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)
from umi import umi
from umi.test_UmiExtractor import exampleForwardRecord, exampleReverseRecord, adapterSequences, topUmi, bottomUmi, targetSequence
from test_conseq import parser, umiArgs, umiFiles
from test_conseq import parsedUmiArgs as args

def test__umi__main(args, topUmi, bottomUmi, exampleForwardRecord, exampleReverseRecord):
    umi.main(args)
    binPath = args["output"] + "bins/"
    assert os.path.isdir(binPath)
    files = os.listdir(binPath)
    assert len(files) == 1
    file = files[0]
    assert file == "targetSequenceBin0.fastq"
    targetSequenceRecordOutputs = list(SeqIO.parse(binPath + file, "fastq"))
    assert len(targetSequenceRecordOutputs) == 2
    if targetSequenceRecordOutputs[0].id == "forward": targetForwardRecordOutput, targetReverseRecordOutput = targetSequenceRecordOutputs
    else: targetReverseRecordOutput, targetForwardRecordOutput = targetSequenceRecordOutputs
    targetForwardRecord = exampleForwardRecord[166:-162]
    targetReverseRecord = exampleReverseRecord[162:-166].reverse_complement()
    assert targetForwardRecordOutput.id == "forward"
    assert targetReverseRecordOutput.id == "reverse"
    assert len(targetForwardRecordOutput) == 200
    assert targetForwardRecordOutput.letter_annotations["phred_quality"] == [40 for _ in range(200)]
    assert len(targetReverseRecordOutput) == 200
    assert targetReverseRecordOutput.letter_annotations["phred_quality"] == [40 for _ in range(200)]

    dataAnalysisPath = args["output"] + "data_analysis/"
    assert os.path.isdir(dataAnalysisPath)
    files = sorted(os.listdir(dataAnalysisPath))
    assert len(files) == 4
    assert os.path.isfile(dataAnalysisPath + "starcode_output_for_top_umis.csv")
    assert os.path.isfile(dataAnalysisPath + "starcode_output_for_bottom_umis.csv")
    assert os.path.isfile(dataAnalysisPath + "chimera_summary_of_starcode_matches.csv")
    assert os.path.isfile(dataAnalysisPath + "read_error_summary.csv")
    chimeraOutput = pd.read_csv(dataAnalysisPath + "chimera_summary_of_starcode_matches.csv")
    assert list(chimeraOutput.columns) == ["top UMI", "bottom UMI", "Number of Reads", "Read Identifiers", "Not Chimera",]
    readErrorOutput = pd.read_csv(dataAnalysisPath + "read_error_summary.csv")
    assert list(readErrorOutput.columns) == ["Read ID","Adapter not found", "Top UMI not found", "Bottom UMI not found", "Target Sequence not found"]

def test__umi__main_fails_when_no_umis_found(args):
    badRecords = [SeqRecord(Seq("A"*200),id=str(i)) for i in range(10)]
    args["input"] = badRecords
    errorOutput = "All provided reads were rejected because no UMIs or target sequences were identified. Please see the 'data_analysis/read_error_summary.csv' file in the output for information on why all reads were rejected."
    with pytest.raises(RuntimeError, match=re.escape(errorOutput)):
        umi.main(args)


def test__umi__starcode():
    originalUmis = [
        "AAAAAAAAA",
        "TTTTTTTTT",
        "CCCCCCCCC",
    ]
    umiText = []
    for _ in range(3):
        for i in range(len(originalUmis)):
            umiText.append(originalUmis[i])

    umiText.append(originalUmis[0])
    umiText.append(originalUmis[0])
    umiText.append(originalUmis[1])

    umiToReadIndicesDict = {
        originalUmis[0]:{1, 4, 7, 10, 11},
        originalUmis[1]: {2, 5, 8, 12},
        originalUmis[2]: {3, 6, 9},
    }
    umiToReadIndicesDictOutput = umi.starcode(umiText)
    assert umiToReadIndicesDictOutput == umiToReadIndicesDict

def test__umi__starcode_does_not_fail_when_only_one_instance_of_umi_found():
    umiText = "AAAAAAAAA"
    umiTextOutput = umi.starcode([umiText])
    assert len(umiTextOutput) == 1
    assert umiTextOutput[umiText] == set([1])



