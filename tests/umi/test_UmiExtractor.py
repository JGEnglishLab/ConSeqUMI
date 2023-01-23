import pytest
import sys
sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
from umi.UmiExtractor import UmiExtractor
from general import genomicFunctions
from unittest.mock import Mock
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
def adapterSeqs():
    adapterSeqs = {
        "topFrontAdapter":"GAGTGTGGCTCTTCGGAT",
        "topBackAdapter":"CACCTTCGTGACTTCCCATT",
        "bottomFrontAdapter":"GTGGGACTGCTGATGACGACTGAT",
        "bottomBackAdapter":"GCGATGCAATTTCCTCATTT",
        "topFrontAdapter_reverseComplement":"ATCCGAAGAGCCACACTC",
        "topBackAdapter_reverseComplement":"AATGGGAAGTCACGAAGGTG",
        "bottomFrontAdapter_reverseComplement":"ATCAGTCGTCATCAGCAGTCCCAC",
        "bottomBackAdapter_reverseComplement":"AAATGAGGAAATTGCATCGC",
    }
    return adapterSeqs

@pytest.fixture
def umiExtractorBasic():
    return UmiExtractor()

def test_umi_extractor_set_umi_pattern_with_empty_umi(umiExtractorBasic):
    umiExtractorBasic.set_umi_pattern("")
    assert umiExtractorBasic.umiPattern == "^$"

def test_umi_extractor_set_umi_pattern_with_nonEmpty_umi(umiExtractorBasic, umiPattern, umiPatternRegEx):
    umiExtractorBasic.set_umi_pattern(umiPattern)
    assert umiExtractorBasic.umiPattern == umiPatternRegEx

def test_umi_extractor_create_linked_adapter(umiExtractorBasic, adapterSeqs):
    testName = "adapterTest"
    linkedAdapter = umiExtractorBasic.create_linked_adapter(
        adapterSeqs["topFrontAdapter"],
        adapterSeqs["topBackAdapter"],
        name=testName,
    )
    assert linkedAdapter.front_adapter.sequence == adapterSeqs["topFrontAdapter"]
    assert linkedAdapter.back_adapter.sequence == adapterSeqs["topBackAdapter"]
    assert linkedAdapter.name == testName

def test_umi_extractor_set_universal_top_and_bottom_linked_adapters(umiExtractorBasic, adapterSeqs):
    umiExtractorBasic.set_universal_top_and_bottom_linked_adapters(
        adapterSeqs["topFrontAdapter"],
        adapterSeqs["topBackAdapter"],
        adapterSeqs["bottomFrontAdapter"],
        adapterSeqs["bottomBackAdapter"],
    )
    assert umiExtractorBasic.topLinkedAdapter.front_adapter.sequence == adapterSeqs["topFrontAdapter"]
    assert umiExtractorBasic.topLinkedAdapter.back_adapter.sequence == adapterSeqs["topBackAdapter"]
    assert umiExtractorBasic.topLinkedAdapter.name == "top"

    assert umiExtractorBasic.bottomLinkedAdapter.front_adapter.sequence == adapterSeqs["bottomFrontAdapter"]
    assert umiExtractorBasic.bottomLinkedAdapter.back_adapter.sequence == adapterSeqs["bottomBackAdapter"]
    assert umiExtractorBasic.bottomLinkedAdapter.name == "bottom"

    assert umiExtractorBasic.topLinkedAdapter_reverseComplement.front_adapter.sequence == adapterSeqs["topBackAdapter_reverseComplement"]
    assert umiExtractorBasic.topLinkedAdapter_reverseComplement.back_adapter.sequence == adapterSeqs["topFrontAdapter_reverseComplement"]
    assert umiExtractorBasic.topLinkedAdapter_reverseComplement.name == "top_reverseComplement"

    assert umiExtractorBasic.bottomLinkedAdapter_reverseComplement.front_adapter.sequence == adapterSeqs["bottomBackAdapter_reverseComplement"]
    assert umiExtractorBasic.bottomLinkedAdapter_reverseComplement.back_adapter.sequence == adapterSeqs["bottomFrontAdapter_reverseComplement"]
    assert umiExtractorBasic.bottomLinkedAdapter_reverseComplement.name == "bottom_reverseComplement"


def test_umi_extractor_set_universal_top_and_bottom_linked_adapters_accounts_for_only_standard_nucleotides(umiExtractorBasic, adapterSeqs):
    errorOutput = "Provided Sequence must only contain standard nucleotides (A, T, C, G)"
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        umiExtractorBasic.set_universal_top_and_bottom_linked_adapters(
            adapterSeqs["topFrontAdapter"] + "R",
            adapterSeqs["topBackAdapter"],
            adapterSeqs["bottomFrontAdapter"],
            adapterSeqs["bottomBackAdapter"],
        )


def test_umi_extractor_initialization_with_umi_pattern(umiPattern, umiPatternRegEx):
    umiExtractor = UmiExtractor(umiPattern=umiPattern)
    assert umiExtractor.umiPattern == umiPatternRegEx

def test_umi_extractor_initialization_with_adapter_sequences(adapterSeqs):
    umiExtractor = UmiExtractor(
        topFrontAdapter=adapterSeqs["topFrontAdapter"],
        topBackAdapter=adapterSeqs["topBackAdapter"],
        bottomFrontAdapter=adapterSeqs["bottomFrontAdapter"],
        bottomBackAdapter=adapterSeqs["bottomBackAdapter"],
    )
    assert umiExtractor.topLinkedAdapter.front_adapter.sequence == adapterSeqs["topFrontAdapter"]
    assert umiExtractor.topLinkedAdapter.back_adapter.sequence == adapterSeqs["topBackAdapter"]
    assert umiExtractor.topLinkedAdapter.name == "top"

    assert umiExtractor.bottomLinkedAdapter.front_adapter.sequence == adapterSeqs["bottomFrontAdapter"]
    assert umiExtractor.bottomLinkedAdapter.back_adapter.sequence == adapterSeqs["bottomBackAdapter"]
    assert umiExtractor.bottomLinkedAdapter.name == "bottom"

    assert umiExtractor.topLinkedAdapter_reverseComplement.front_adapter.sequence == adapterSeqs["topBackAdapter_reverseComplement"]
    assert umiExtractor.topLinkedAdapter_reverseComplement.back_adapter.sequence == adapterSeqs["topFrontAdapter_reverseComplement"]
    assert umiExtractor.topLinkedAdapter_reverseComplement.name == "top_reverseComplement"

    assert umiExtractor.bottomLinkedAdapter_reverseComplement.front_adapter.sequence == adapterSeqs["bottomBackAdapter_reverseComplement"]
    assert umiExtractor.bottomLinkedAdapter_reverseComplement.back_adapter.sequence == adapterSeqs["bottomFrontAdapter_reverseComplement"]
    assert umiExtractor.bottomLinkedAdapter_reverseComplement.name == "bottom_reverseComplement"

@pytest.fixture
def umiExtractor(umiPattern, adapterSeqs):
    umiExtractor = UmiExtractor(
        umiPattern=umiPattern,
        topFrontAdapter=adapterSeqs["topFrontAdapter"],
        topBackAdapter=adapterSeqs["topBackAdapter"],
        bottomFrontAdapter=adapterSeqs["bottomFrontAdapter"],
        bottomBackAdapter=adapterSeqs["bottomBackAdapter"],
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

@pytest.fixture
def exampleForwardRecord(adapterSeqs, topUmi, bottomUmi, targetSequence):
    fillerSeq1 = "A" * 110
    fillerSeq2 = "T" * 100
    exampleForwardSequence = "".join(
        [
            fillerSeq1,
            adapterSeqs["topFrontAdapter"],
            topUmi,
            adapterSeqs["topBackAdapter"],
            targetSequence,
            adapterSeqs["bottomBackAdapter_reverseComplement"],
            genomicFunctions.find_reverse_complement(bottomUmi),
            adapterSeqs["bottomFrontAdapter_reverseComplement"],
            fillerSeq2,
        ]
    )
    return SeqRecord(Seq(exampleForwardSequence), id="forward")

@pytest.fixture
def exampleReverseRecord(exampleForwardRecord):
    return SeqRecord(exampleForwardRecord.seq.reverse_complement, id="reverse")

def test_umi_extractor_find_matches_of_adapters_in_sequence(umiExtractor, exampleForwardRecord):
    sequence = str(exampleForwardRecord.seq)
    topMatch, bottomMatch = umiExtractor.find_matches_of_adapters_in_sequence(sequence)
    assert topMatch is not None
    assert bottomMatch is not None

def test_umi_extractor_find_matches_of_adapters_in_sequence_when_no_match_found(umiExtractor, exampleForwardRecord):
    sequence = str(exampleForwardRecord.seq)
    topErrorSequence = "A"*200 + sequence[200:]
    topMatchError, bottomMatch = umiExtractor.find_matches_of_adapters_in_sequence(topErrorSequence)
    assert topMatchError is None
    assert bottomMatch is not None

    bottomErrorSequence = sequence[:-200] + "A"*200
    topMatch, bottomMatchError = umiExtractor.find_matches_of_adapters_in_sequence(bottomErrorSequence)
    assert topMatch is not None
    assert bottomMatchError is None


#def test_umi_extractor_extract_umis_and_target_sequence_from_record_of_forward_sequence(umiExtractor, exampleForwardRecord, topUmi, bottomUmi, targetSequence):
#    topUmiOutput, bottomUmiOutput, targetSequenceOutput = umiExtractor.extract_umis_and_target_sequence_from_record(exampleForwardRecord)
