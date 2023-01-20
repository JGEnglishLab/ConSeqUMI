import pytest
import sys
sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
from umi.UmiExtractor import UmiExtractor
from unittest.mock import Mock

def test_umi_extractor_initialization():
    ue = UmiExtractor()
    assert ue.umiPattern == "^$"

@pytest.fixture
def adapterSeqs():
    adapterSeqs = {
        "forwardAdapterFront":"CAAGCAGAAGACGGCATACGAGAT",
        "forwardAdapterBack":"AGRGTTYGATYMTGGCTCAG",
        "reverseAdapterFront":"AATGATACGGCGACCACCGAGATC",
        "reverseAdapterBack":"CGACATCGAGGTGCCAAAC",
        "forwardAdapterFront_reverseComplement":"ATCTCGTATGCCGTCTTCTGCTTG",
        "forwardAdapterBack_reverseComplement":"CTGAGCCAKRATCRAACYCT",
        "reverseAdapterFront_reverseComplement":"GATCTCGGTGGTCGCCGTATCATT",
        "reverseAdapterBack_reverseComplement":"GTTTGGCACCTCGATGTCG",
    }
    return adapterSeqs

@pytest.fixture
def ue():
    ue = UmiExtractor()
    ue.umiPattern = Mock(return_value="^[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{6}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}$")
    return ue

def test_umi_extractor_create_linked_adapter(ue, adapterSeqs):
    testName = "adapterTest"
    linkedAdapter = ue.create_linked_adapter(
        adapterSeqs["forwardAdapterFront"],
        adapterSeqs["forwardAdapterBack"],
        name=testName,
    )
    assert linkedAdapter.front_adapter.sequence == adapterSeqs["forwardAdapterFront"]
    assert linkedAdapter.back_adapter.sequence == adapterSeqs["forwardAdapterBack"]
    assert linkedAdapter.name == testName
