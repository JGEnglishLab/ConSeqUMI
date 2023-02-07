import pytest
import random
import os
import pandas as pd
from Bio import SeqIO

import test__umi
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



