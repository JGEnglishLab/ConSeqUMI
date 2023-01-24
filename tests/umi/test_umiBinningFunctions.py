import pytest
import sys
sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
from umi import umiBinningFunctions
import pandas as pd
#from unittest.mock import Mock
#import re
#import random
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord

@pytest.fixture
def topUmiToReadIndices():
    topUmiToReadIndices = {
        "AAAATTTT":{1,2,3,4,5,6,7,8,9,10},
        "AATTAATT":{11,12,13,14,15},
        "TTTTAAAA":{16,17,18,19,20},
        "TTAATTAA":{21},
    }
    return topUmiToReadIndices

@pytest.fixture
def bottomUmiToReadIndices():
    bottomUmiToReadIndices = {
        "CCCCGGGG":{1,2,3,4,5,6,11,12},
        "CCGGCCGG":{7,8,9,10,13,14,15,21},
        "GGGGCCCC":{16,17,18,19,20},
    }
    return bottomUmiToReadIndices

@pytest.fixture
def pairedUmiToReadIndices():
    pairedUmiToReadIndices = {
        "AAAATTTTCCCCGGGG":{1,2,3,4,5,6,7,8,9,10},
        "TTTTAAAAGGGGCCCC":{16,17,18,19,20},
        "AATTAATTCCGGCCGG":{13,14,15},
    }
    return pairedUmiToReadIndices

@pytest.fixture
def allUmiPairStructure():
    topUmis = [
        "AAAATTTT",
        "TTTTAAAA",
        "AAAATTTT",
        "AATTAATT",
        "AATTAATT",
        "TTAATTAA",
    ]
    bottomUmis = [
        "CCCCGGGG",
        "GGGGCCCC",
        "CCGGCCGG",
        "CCGGCCGG",
        "CCCCGGGG",
        "CCGGCCGG",
    ]
    readIndices = [
        {1,2,3,4,5,6},
        {16,17,18,19,20},
        {7,8,9,10},
        {13,14,15},
        {11,12},
        {21},
    ]
    return [topUmis, bottomUmis, readIndices]

def test_umi_binning_functions_pair_top_and_bottom_umi_by_matching_reads(topUmiToReadIndices, bottomUmiToReadIndices, allUmiPairStructure):
    topUmis, bottomUmis, readIndices = allUmiPairStructure
    topUmisOutput, bottomUmisOutput, readIndicesOutput = umiBinningFunctions.pair_top_and_bottom_umi_by_matching_reads(topUmiToReadIndices, bottomUmiToReadIndices)
    assert topUmis == topUmisOutput
    assert bottomUmis == bottomUmisOutput
    assert readIndices == readIndicesOutput
