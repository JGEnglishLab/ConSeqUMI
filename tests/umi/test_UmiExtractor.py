import pytest
import sys
sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
from umi.UmiExtractor import UmiExtractor
from unittest.mock import Mock
import re

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
    with pytest.raises(TypeError, match=re.escape(errorOutput)):
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
