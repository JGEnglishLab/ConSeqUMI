import pytest
from umi import umiBinningFunctions
import pandas as pd

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

@pytest.fixture
def chimeraIndices():
    return [2, 4, 5]

def test_umi_binning_functions_identify_chimera_indices(allUmiPairStructure, chimeraIndices):
    topUmis, bottomUmis, readIndices = allUmiPairStructure
    chimeraIndicesOutput = umiBinningFunctions.identify_chimera_indices(topUmis, bottomUmis)
    assert chimeraIndicesOutput == chimeraIndices

@pytest.fixture
def pairedUmiToReadIndices():
    pairedUmiToReadIndices = {
        "AAAATTTTCCCCGGGG":{1,2,3,4,5,6},
        "TTTTAAAAGGGGCCCC":{16,17,18,19,20},
        "AATTAATTCCGGCCGG":{13,14,15},
    }
    return pairedUmiToReadIndices

def test_umi_binning_functions_remove_chimeras_from_umi_pairs_and_return_paired_umi_to_read_indices_dict(allUmiPairStructure, chimeraIndices, pairedUmiToReadIndices):
    pairedUmiToReadIndicesOutput = umiBinningFunctions.remove_chimeras_from_umi_pairs_and_return_paired_umi_to_read_indices_dict(*allUmiPairStructure, chimeraIndices)
    assert pairedUmiToReadIndicesOutput == pairedUmiToReadIndices

def test_umi_binning_functions_compile_chimera_data_analysis_data_frame(allUmiPairStructure, chimeraIndices):
    topUmis, bottomUmis, readIndices = allUmiPairStructure
    readIndicesString = ["/".join([str(index) for index in sorted(readIndexSet)]) for readIndexSet in readIndices]
    readIndicesLength = [len(readIndexSet) for readIndexSet in readIndices]
    nonChimeraIndicator = [1 if i in chimeraIndices else 0 for i in range(len(topUmis))]
    chimeraDataValues = [
        topUmis,
        bottomUmis,
        readIndicesLength,
        readIndicesString,
        nonChimeraIndicator,
    ]
    chimeraDataColumnNames = [
        "top UMI",
        "bottom UMI",
        "Number of Reads",
        "Read Identifiers",
        "Not Chimera",
    ]
    chimeraData = pd.DataFrame(list(zip(*chimeraDataValues)), columns=chimeraDataColumnNames)
    chimeraDataOutput = umiBinningFunctions.compile_chimera_data_analysis_data_frame(*allUmiPairStructure, chimeraIndices)
    assert chimeraDataOutput.equals(chimeraData)
