import pytest
import importlib
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)
from umi.UmiExtractor import UmiExtractor
from umi import umiExtractionFunctions
import re
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

@pytest.fixture
def umiPattern():
    return "NNNYRNNNYRNNNYRNNNNNNYRNNNYRNNNYRNNN"
@pytest.fixture
def umiPatternRegEx():
    return "^[ACGT][ACGT][ACGT][CT][AG][ACGT][ACGT][ACGT][CT][AG][ACGT][ACGT][ACGT][CT][AG][ACGT][ACGT][ACGT][ACGT][ACGT][ACGT][CT][AG][ACGT][ACGT][ACGT][CT][AG][ACGT][ACGT][ACGT][CT][AG][ACGT][ACGT][ACGT]$"

@pytest.fixture
def adapterSequences():
    adapterSequences = {
        "topFrontAdapter":"GAGTGTGGCTCTTCGGAT",
        "topBackAdapter":"CACCTTCGTGACTTCCCATT",
        "bottomFrontAdapter":"GTGGGACTGCTGATGACGACTGAT",
        "bottomBackAdapter":"GCGATGCAATTTCCTCATTT",
        "topFrontAdapter_reverseComplement":"ATCCGAAGAGCCACACTC",
        "topBackAdapter_reverseComplement":"AATGGGAAGTCACGAAGGTG",
        "bottomFrontAdapter_reverseComplement":"ATCAGTCGTCATCAGCAGTCCCAC",
        "bottomBackAdapter_reverseComplement":"AAATGAGGAAATTGCATCGC",
    }
    return adapterSequences

@pytest.fixture
def umiExtractorBasic():
    return UmiExtractor()

def test__umi_extractor__set_umi_pattern_with_empty_umi(umiExtractorBasic):
    umiExtractorBasic.set_umi_pattern("")
    assert umiExtractorBasic.umiPattern == "^$"

def test__umi_extractor__set_umi_pattern_with_nonEmpty_umi(umiExtractorBasic, umiPattern, umiPatternRegEx):
    umiExtractorBasic.set_umi_pattern(umiPattern)
    assert umiExtractorBasic.umiPattern == umiPatternRegEx

def test__umi_extractor__create_linked_adapter(umiExtractorBasic, adapterSequences):
    testName = "adapterTest"
    linkedAdapter = umiExtractorBasic.create_linked_adapter(
        adapterSequences["topFrontAdapter"],
        adapterSequences["topBackAdapter"],
        name=testName,
    )
    assert linkedAdapter.front_adapter.sequence == adapterSequences["topFrontAdapter"]
    assert linkedAdapter.back_adapter.sequence == adapterSequences["topBackAdapter"]
    assert linkedAdapter.name == testName

def test__umi_extractor__set_universal_top_and_bottom_linked_adapters(umiExtractorBasic, adapterSequences):
    umiExtractorBasic.set_universal_top_and_bottom_linked_adapters(
        adapterSequences["topFrontAdapter"],
        adapterSequences["topBackAdapter"],
        adapterSequences["bottomFrontAdapter"],
        adapterSequences["bottomBackAdapter"],
    )
    assert umiExtractorBasic.topAdapter.front_adapter.sequence == adapterSequences["topFrontAdapter"]
    assert umiExtractorBasic.topAdapter.back_adapter.sequence == adapterSequences["topBackAdapter"]
    assert umiExtractorBasic.topAdapter.name == "top"

    assert umiExtractorBasic.bottomAdapter.front_adapter.sequence == adapterSequences["bottomFrontAdapter"]
    assert umiExtractorBasic.bottomAdapter.back_adapter.sequence == adapterSequences["bottomBackAdapter"]
    assert umiExtractorBasic.bottomAdapter.name == "bottom"

def test__umi_extractor__set_universal_top_and_bottom_linked_adapters_accounts_for_only_standard_nucleotides(umiExtractorBasic, adapterSequences):
    errorOutput = "Provided Adapter Sequences must only contain standard nucleotides (A, T, C, G)"
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        umiExtractorBasic.set_universal_top_and_bottom_linked_adapters(
            adapterSequences["topFrontAdapter"] + "R",
            adapterSequences["topBackAdapter"],
            adapterSequences["bottomFrontAdapter"],
            adapterSequences["bottomBackAdapter"],
        )


def test__umi_extractor__initialization_with_umi_pattern(umiPattern, umiPatternRegEx):
    umiExtractor = UmiExtractor(umiPattern=umiPattern)
    assert umiExtractor.umiPattern == umiPatternRegEx

def test__umi_extractor__initialization_with_adapter_sequences(adapterSequences):
    umiExtractor = UmiExtractor(
        topFrontAdapter=adapterSequences["topFrontAdapter"],
        topBackAdapter=adapterSequences["topBackAdapter"],
        bottomFrontAdapter=adapterSequences["bottomFrontAdapter"],
        bottomBackAdapter=adapterSequences["bottomBackAdapter"],
    )
    assert umiExtractor.topAdapter.front_adapter.sequence == adapterSequences["topFrontAdapter"]
    assert umiExtractor.topAdapter.back_adapter.sequence == adapterSequences["topBackAdapter"]
    assert umiExtractor.topAdapter.name == "top"

    assert umiExtractor.bottomAdapter.front_adapter.sequence == adapterSequences["bottomFrontAdapter"]
    assert umiExtractor.bottomAdapter.back_adapter.sequence == adapterSequences["bottomBackAdapter"]
    assert umiExtractor.bottomAdapter.name == "bottom"

@pytest.fixture
def umiExtractor(umiPattern, adapterSequences):
    umiExtractor = UmiExtractor(
        umiPattern=umiPattern,
        topFrontAdapter=adapterSequences["topFrontAdapter"],
        topBackAdapter=adapterSequences["topBackAdapter"],
        bottomFrontAdapter=adapterSequences["bottomFrontAdapter"],
        bottomBackAdapter=adapterSequences["bottomBackAdapter"],
    )
    return umiExtractor

@pytest.fixture
def targetSequence():
    random.seed(0)
    targetSequence = "".join(random.choices("ATGC", k=200))
    return targetSequence

@pytest.fixture
def topUmi(): return "ATATGTACTATTGTATAC"

@pytest.fixture
def bottomUmi(): return "CAATGTGTCGGCATAGGG"

def make_string_to_fastq_record(string, id):
    sequence = Seq(string)
    id = str(id)
    phred_quality = [40 for j in range(len(string))]
    return SeqRecord(sequence, id=id, letter_annotations={"phred_quality":phred_quality})

@pytest.fixture
def exampleForwardRecord(adapterSequences, topUmi, bottomUmi, targetSequence):
    fillerSeq1 = "A" * 110
    fillerSeq2 = "T" * 100
    exampleForwardSequence = "".join(
        [
            fillerSeq1,
            adapterSequences["topFrontAdapter"],
            topUmi,
            adapterSequences["topBackAdapter"],
            targetSequence,
            adapterSequences["bottomBackAdapter_reverseComplement"],
            umiExtractionFunctions.find_reverse_complement(bottomUmi),
            adapterSequences["bottomFrontAdapter_reverseComplement"],
            fillerSeq2,
        ]
    )
    return make_string_to_fastq_record(exampleForwardSequence, id="forward")

@pytest.fixture
def exampleReverseRecord(exampleForwardRecord):
    return make_string_to_fastq_record(exampleForwardRecord.seq.reverse_complement(), id="reverse")

def test__umi_extractor__find_matches_of_adapters_in_sequence(umiExtractor, exampleForwardRecord):
    sequence = str(exampleForwardRecord.seq)
    topMatch, bottomMatch = umiExtractor.find_matches_of_adapters_in_sequence(sequence)
    assert topMatch is not None
    assert bottomMatch is not None

def test__umi_extractor__find_matches_of_adapters_in_sequence_when_no_match_found(umiExtractor, exampleForwardRecord):
    sequence = str(exampleForwardRecord.seq)
    topErrorSequence = "A"*200 + sequence[200:]
    topMatchError, bottomMatch = umiExtractor.find_matches_of_adapters_in_sequence(topErrorSequence)
    assert topMatchError is None
    assert bottomMatch is not None

    bottomErrorSequence = sequence[:-200] + "A"*200
    topMatch, bottomMatchError = umiExtractor.find_matches_of_adapters_in_sequence(bottomErrorSequence)
    assert topMatch is not None
    assert bottomMatchError is None

def test__umi_extractor__extract_umis_and_target_sequence_from_read__record_of_forward_sequence(umiExtractor, exampleForwardRecord, topUmi, bottomUmi, targetSequence):
    topUmiOutput, bottomUmiOutput, targetSequenceRecordOutput = umiExtractor.extract_umis_and_target_sequence_from_read(exampleForwardRecord)
    assert topUmiOutput == topUmi
    assert bottomUmiOutput == bottomUmi
    assert str(targetSequenceRecordOutput.seq) == targetSequence
    assert targetSequenceRecordOutput.id == exampleForwardRecord.id
    assert targetSequenceRecordOutput.letter_annotations == exampleForwardRecord[:len(targetSequence)].letter_annotations

def test__umi_extractor__extract_umis_and_target_sequence_from_read__record_of_reverse_sequence(umiExtractor, exampleReverseRecord, topUmi, bottomUmi, targetSequence):
    topUmiOutput, bottomUmiOutput, targetSequenceRecordOutput = umiExtractor.extract_umis_and_target_sequence_from_read(exampleReverseRecord)
    assert topUmiOutput == topUmi
    assert bottomUmiOutput == bottomUmi
    assert str(targetSequenceRecordOutput.seq) == targetSequence
    assert targetSequenceRecordOutput.id == exampleReverseRecord.id
    assert targetSequenceRecordOutput.letter_annotations == exampleReverseRecord[:len(targetSequence)].letter_annotations

@pytest.fixture
def exampleForwardRecord_withTopUmiNotFound(exampleForwardRecord):
    exampleForwardSequence = str(exampleForwardRecord.seq)
    exampleForwardSequence_withTopUmiNotFound = "A"*200 + exampleForwardSequence[200:]
    exampleForwardRecord_withTopUmiNotFound = SeqRecord(Seq(exampleForwardSequence_withTopUmiNotFound), id="forward top error")
    return exampleForwardRecord_withTopUmiNotFound

@pytest.fixture
def exampleForwardRecord_withBottomUmiNotFound(exampleForwardRecord):
    exampleForwardSequence = str(exampleForwardRecord.seq)
    exampleForwardSequence_withBottomUmiNotFound = exampleForwardSequence[:-200] + "A"*200
    exampleForwardRecord_withBottomUmiNotFound = SeqRecord(Seq(exampleForwardSequence_withBottomUmiNotFound), id="forward Bottom error")
    return exampleForwardRecord_withBottomUmiNotFound

def test__umi_extractor__extract_umis_and_target_sequence_from_record_errors_when_no_umi_found_in_top_or_bottom(umiExtractor, exampleForwardRecord_withTopUmiNotFound, exampleForwardRecord_withBottomUmiNotFound):
    topUmiOutput, bottomUmiOutput, targetSequenceRecordOutput = umiExtractor.extract_umis_and_target_sequence_from_read(exampleForwardRecord_withTopUmiNotFound)
    assert topUmiOutput == ""
    assert bottomUmiOutput == ""
    assert str(targetSequenceRecordOutput.seq) == ""
    assert targetSequenceRecordOutput.name == "adapter not found"

    topUmiOutput, bottomUmiOutput, targetSequenceRecordOutput = umiExtractor.extract_umis_and_target_sequence_from_read(exampleForwardRecord_withBottomUmiNotFound)
    assert topUmiOutput == ""
    assert bottomUmiOutput == ""
    assert str(targetSequenceRecordOutput.seq) == ""
    assert targetSequenceRecordOutput.name == "adapter not found"

def test__umi_extractor__extract_umis_and_target_sequences_from_all_records(umiExtractor, exampleForwardRecord, exampleReverseRecord, topUmi, bottomUmi, targetSequence):
    topUmisOutput, bottomUmisOutput, targetSequenceRecordsOutput = umiExtractor.extract_umis_and_target_sequences_from_all_records([exampleForwardRecord, exampleReverseRecord])
    assert set(topUmisOutput) == {topUmi}
    assert set(bottomUmisOutput) == {bottomUmi}
    assert str(targetSequenceRecordsOutput[0].seq) == targetSequence
    assert str(targetSequenceRecordsOutput[1].seq) == targetSequence
    assert targetSequenceRecordsOutput[0].id == exampleForwardRecord.id
    assert targetSequenceRecordsOutput[1].id == exampleReverseRecord.id
