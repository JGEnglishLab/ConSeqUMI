import pytest
from umi import umiBinningFunctions
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from unittest.mock import Mock

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

def test__umi_binning_functions__pair_top_and_bottom_umi_by_matching_reads(topUmiToReadIndices, bottomUmiToReadIndices, allUmiPairStructure):
    topUmis, bottomUmis, readIndices = allUmiPairStructure
    topUmisOutput, bottomUmisOutput, readIndicesOutput = umiBinningFunctions.pair_top_and_bottom_umi_by_matching_reads(topUmiToReadIndices, bottomUmiToReadIndices)
    assert topUmis == topUmisOutput
    assert bottomUmis == bottomUmisOutput
    assert readIndices == readIndicesOutput

@pytest.fixture
def chimeraIndices():
    return [2, 4, 5]

@pytest.fixture
def umiExtractionCorrectOutput():
    return ["AAAATTTT"], ["CCCCGGGG"], [SeqRecord(Seq("AAA"),id="correctOutput")]

def test__umi_binning_functions__identify_indices_of_reads_missing_key_values__identifies_when_adapter_not_found():
    pass

def test__umi_binning_functions__identify_chimera_indices(allUmiPairStructure, chimeraIndices):
    topUmis, bottomUmis, readIndices = allUmiPairStructure
    chimeraIndicesOutput = umiBinningFunctions.identify_chimera_indices(topUmis, bottomUmis)
    assert chimeraIndicesOutput == chimeraIndices

def test__umi_binning_functions__remove_chimeras_from_umi_pairs_and_return_paired_umi_to_read_records_dict(allUmiPairStructure, chimeraIndices):
    targetSequence = "AAA"
    phred_quality = [40 for j in range(len(targetSequence))]
    recordsInput = [SeqRecord(Seq(targetSequence),id=str(i), letter_annotations={"phred_quality":phred_quality}) for i in range(1,21)]
    recordsOutput = recordsInput[:]
    binNumbers = [
        [1,2,3,4,5,6],
        [16,17,18,19,20],
        [13,14,15],
    ]
    for i in binNumbers[0]:
        recordsOutput[i-1].description = f"Top UMI: AAAATTTT; Bottom UMI: CCCCGGGG; read number: {i}"
    for i in binNumbers[1]:
        recordsOutput[i-1].description = f"Top UMI: TTTTAAAA; Bottom UMI: GGGGCCCC; read number: {i}"
    for i in binNumbers[2]:
        recordsOutput[i-1].description = f"Top UMI: AATTAATT, Bottom UMI: CCGGCCGG; read number: {i}"

    pairedUmiToReadRecords = {
        ("AAAATTTT","CCCCGGGG"):[recordsOutput[i-1] for i in binNumbers[0]],
        ("TTTTAAAA","GGGGCCCC"):[recordsOutput[i-1] for i in binNumbers[1]],
        ("AATTAATT","CCGGCCGG"):[recordsOutput[i-1] for i in binNumbers[2]],
    }
    pairedUmiToReadRecordsOutput = umiBinningFunctions.remove_chimeras_from_umi_pairs_and_return_paired_umi_to_read_records_dict(*allUmiPairStructure, chimeraIndices, recordsInput)
    assert pairedUmiToReadRecordsOutput == pairedUmiToReadRecords

def test__umi_binning_functions__compile_chimera_data_analysis_data_frame(allUmiPairStructure, chimeraIndices):
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
