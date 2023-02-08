import pytest
import random
import os
import pandas as pd
from Bio import SeqIO

import umi.test__umi
from umi import umi
umi.SCOMMAND = ['starcode', '--seq-id', '-q']
from umi.test_UmiExtractor import exampleForwardRecord, exampleReverseRecord, adapterSequences, topUmi, bottomUmi, targetSequence
from test_conseq import parser, umiArgs, files
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

#os.listdir(
#os.path.isfile(name)



