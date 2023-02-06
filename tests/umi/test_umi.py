import pytest
import random
import os
import pandas as pd
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
    #assert len(files) == 1

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



