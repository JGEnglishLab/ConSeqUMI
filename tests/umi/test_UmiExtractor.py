import pytest
import sys
sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
from umi.UmiExtractor import UmiExtractor
from unittest.mock import Mock

def test_umi_extractor_initialization():
    umiExtractor = UmiExtractor()
    assert umiExtractor.umiPattern == "^$"

@pytest.fixture
def adapterSeqs():
    adapterSeqs = {
        "topAdapterFront":"CAAGCAGAAGACGGCATACGAGAT",
        "topAdapterBack":"AGRGTTYGATYMTGGCTCAG",
        "bottomAdapterFront":"AATGATACGGCGACCACCGAGATC",
        "bottomAdapterBack":"CGACATCGAGGTGCCAAAC",
        "topAdapterFront_reverseComplement":"ATCTCGTATGCCGTCTTCTGCTTG",
        "topAdapterBack_reverseComplement":"CTGAGCCAKRATCRAACYCT",
        "bottomAdapterFront_reverseComplement":"GATCTCGGTGGTCGCCGTATCATT",
        "bottomAdapterBack_reverseComplement":"GTTTGGCACCTCGATGTCG",
    }
    return adapterSeqs

@pytest.fixture
def umiExtractorBasic():
    umiExtractorBasic = UmiExtractor()
    umiExtractorBasic.umiPattern = Mock(return_value="^[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{6}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}$")
    return umiExtractorBasic

def test_umi_extractor_create_linked_adapter(umiExtractorBasic, adapterSeqs):
    testName = "adapterTest"
    linkedAdapter = umiExtractorBasic.create_linked_adapter(
        adapterSeqs["topAdapterFront"],
        adapterSeqs["topAdapterBack"],
        name=testName,
    )
    assert linkedAdapter.front_adapter.sequence == adapterSeqs["topAdapterFront"]
    assert linkedAdapter.back_adapter.sequence == adapterSeqs["topAdapterBack"]
    assert linkedAdapter.name == testName

def test_umi_extractor_set_universal_top_and_bottom_linked_adapters(umiExtractorBasic, adapterSeqs):
    umiExtractorBasic.set_universal_top_and_bottom_linked_adapters(
        adapterSeqs["topAdapterFront"],
        adapterSeqs["topAdapterBack"],
        adapterSeqs["bottomAdapterFront"],
        adapterSeqs["bottomAdapterBack"],
    )
    assert umiExtractorBasic.topLinkedAdapter.front_adapter.sequence == adapterSeqs["topAdapterFront"]
    assert umiExtractorBasic.topLinkedAdapter.back_adapter.sequence == adapterSeqs["topAdapterBack"]
    assert umiExtractorBasic.topLinkedAdapter.name == "top"

    assert umiExtractorBasic.bottomLinkedAdapter.front_adapter.sequence == adapterSeqs["bottomAdapterFront"]
    assert umiExtractorBasic.bottomLinkedAdapter.back_adapter.sequence == adapterSeqs["bottomAdapterBack"]
    assert umiExtractorBasic.bottomLinkedAdapter.name == "bottom"

    assert umiExtractorBasic.topLinkedAdapter_reverseComplement.front_adapter.sequence == adapterSeqs["topAdapterBack_reverseComplement"]
    assert umiExtractorBasic.topLinkedAdapter_reverseComplement.back_adapter.sequence == adapterSeqs["topAdapterFront_reverseComplement"]
    assert umiExtractorBasic.topLinkedAdapter_reverseComplement.name == "top_reverseComplement"

    assert umiExtractorBasic.bottomLinkedAdapter_reverseComplement.front_adapter.sequence == adapterSeqs["bottomAdapterBack_reverseComplement"]
    assert umiExtractorBasic.bottomLinkedAdapter_reverseComplement.back_adapter.sequence == adapterSeqs["bottomAdapterFront_reverseComplement"]
    assert umiExtractorBasic.bottomLinkedAdapter_reverseComplement.name == "bottom_reverseComplement"
