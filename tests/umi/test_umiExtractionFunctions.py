import pytest
from cutadapt.parser import FrontAdapter, BackAdapter, LinkedAdapter

import sys
import os

srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src/ConSeqUMI"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)
from umi import umiExtractionFunctions
import re


@pytest.fixture
def allIUPAC_nucleotides():
    return "ACGTWSRYKMBDHVN"


def test__umi_extraction_functions__convert_IUPAC_to_regular_expression(
    allIUPAC_nucleotides,
):
    nucleotidesInRegExFormat = "ACGT[AT][CG][AG][CT][GT][AC][CGT][AGT][ACT][ACG][ACGT]"
    assert (
        umiExtractionFunctions.convert_IUPAC_to_regular_expression(allIUPAC_nucleotides)
        == nucleotidesInRegExFormat
    )


def test__umi_extraction_functions__is_IUPAC_nucleotide(allIUPAC_nucleotides):
    assert umiExtractionFunctions.is_IUPAC_nucleotide(allIUPAC_nucleotides)


def test__umi_extraction_functions__is_IUPAC_nucleotide_rejects_inappropriate_character(
    allIUPAC_nucleotides,
):
    allIUPAC_nucleotides_withError = allIUPAC_nucleotides + "@"
    assert not umiExtractionFunctions.is_IUPAC_nucleotide(
        allIUPAC_nucleotides_withError
    )


def test__umi_extraction_functions__find_reverse_complement(allIUPAC_nucleotides):
    allIUPAC_nucleotides_reverseComplement = "NBDHVKMRYSWACGT"
    assert (
        umiExtractionFunctions.find_reverse_complement(allIUPAC_nucleotides)
        == allIUPAC_nucleotides_reverseComplement
    )


def test__umi_extraction_functions__find_reverse_complement_non_nucleotide_character_error():
    nonNucleotideCharacter = "@"
    errorOutput = "Provided sequence contains non-nucleotide characters"
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        umiExtractionFunctions.find_reverse_complement(nonNucleotideCharacter)


def test__umi_extraction_functions__find_reverse_complement_works_regardless_of_case():
    standardNucleotides_lowerCase = "atcg"
    standardNucleotides_reverseComplement = "CGAT"
    assert (
        umiExtractionFunctions.find_reverse_complement(standardNucleotides_lowerCase)
        == standardNucleotides_reverseComplement
    )


def test__umi_extraction_functions__extract_top_and_bottom_of_sequence():
    seqExtractionLength = 200
    topSeq = "A" * seqExtractionLength
    middleSeq = "T"
    bottomSeq = "C" * seqExtractionLength
    bottomSeq_reverseComplement = umiExtractionFunctions.find_reverse_complement(
        bottomSeq
    )
    sequence = topSeq + middleSeq + bottomSeq
    (
        topSeqOutput,
        bottomSeqOutput,
    ) = umiExtractionFunctions.extract_top_and_bottom_of_sequence(sequence)
    assert topSeqOutput == topSeq
    assert bottomSeqOutput == bottomSeq_reverseComplement


def test__umi_extraction_functions__find_index_at_end_of_back_adapter__with_front_adapter():
    frontAdapterSequence = "ATCGATCG"
    backAdapterSequence = "AAAATTTT"
    umi = "C" * 20
    sequenceWithFrontAdapter = (
        "G" * 10 + frontAdapterSequence + umi + backAdapterSequence + "G" * 10
    )
    frontAdapter = FrontAdapter(frontAdapterSequence, max_errors=0.2, min_overlap=11)
    backAdapter = BackAdapter(backAdapterSequence, max_errors=0.2, min_overlap=11)
    linkedAdapter = LinkedAdapter(
        frontAdapter,
        backAdapter,
        name="test",
        front_required=True,
        back_required=True,
    )
    match = linkedAdapter.match_to(sequenceWithFrontAdapter)
    indexAtEndOfBackAdapter = len(sequenceWithFrontAdapter) - 10
    indexAtEndOfBackAdapterOutput = (
        umiExtractionFunctions.find_index_at_end_of_back_adapter(match)
    )
    assert indexAtEndOfBackAdapter == indexAtEndOfBackAdapterOutput


def test__umi_extraction_functions__find_index_at_end_of_back_adapter__with_front_adapter():
    frontAdapterSequence = "ATCGATCG"
    backAdapterSequence = "AAAATTTT"
    umi = "C" * 20
    sequenceWithoutFrontAdapter = "G" * 10 + umi + backAdapterSequence + "G" * 10
    frontAdapter = FrontAdapter(frontAdapterSequence, max_errors=0.2, min_overlap=11)
    backAdapter = BackAdapter(backAdapterSequence, max_errors=0.2, min_overlap=11)
    linkedAdapter = LinkedAdapter(
        frontAdapter,
        backAdapter,
        name="test",
        front_required=False,
        back_required=True,
    )
    match = linkedAdapter.match_to(sequenceWithoutFrontAdapter)

    indexAtEndOfBackAdapter = len(sequenceWithoutFrontAdapter) - 10
    indexAtEndOfBackAdapterOutput = (
        umiExtractionFunctions.find_index_at_end_of_back_adapter(match)
    )
    assert indexAtEndOfBackAdapter == indexAtEndOfBackAdapterOutput
