import pytest
from umi import umiExtractionFunctions
import re

@pytest.fixture
def allIUPAC_nucleotides(): return "ACGTWSRYKMBDHVN"


def test_genomic_function_convert_IUPAC_to_regular_expression(allIUPAC_nucleotides):
    nucleotidesInRegExFormat = "ACGT[AT][CG][AG][CT][GT][AC][CGT][AGT][ACT][ACG][ACGT]"
    assert umiExtractionFunctions.convert_IUPAC_to_regular_expression(allIUPAC_nucleotides) == nucleotidesInRegExFormat

def test_genomic_function_is_only_standard_nucleotide():
    standardNucleotides = "ATCG"
    nonStandardNucleotides = "ATCGR"
    assert umiExtractionFunctions.is_only_standard_nucleotide(standardNucleotides)
    assert not umiExtractionFunctions.is_only_standard_nucleotide(nonStandardNucleotides)

def test_genomic_function_is_IUPAC_nucleotide(allIUPAC_nucleotides):
    assert umiExtractionFunctions.is_IUPAC_nucleotide(allIUPAC_nucleotides)

def test_genomic_function_is_IUPAC_nucleotide_rejects_inappropriate_character(allIUPAC_nucleotides):
    allIUPAC_nucleotides_withError = allIUPAC_nucleotides + "@"
    assert not umiExtractionFunctions.is_IUPAC_nucleotide(allIUPAC_nucleotides_withError)

def test_genomic_function_find_reverse_complement(allIUPAC_nucleotides):
    allIUPAC_nucleotides_reverseComplement = "NBDHVKMRYSWACGT"
    assert umiExtractionFunctions.find_reverse_complement(allIUPAC_nucleotides) == allIUPAC_nucleotides_reverseComplement

def test_genomic_function_find_reverse_complement_non_nucleotide_character_error():
    nonNucleotideCharacter = "@"
    errorOutput = "Provided sequence contains non-nucleotide characters"
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        umiExtractionFunctions.find_reverse_complement(nonNucleotideCharacter)

def test_genomic_function_find_reverse_complement_works_regardless_of_case():
    standardNucleotides_lowerCase = "atcg"
    standardNucleotides_reverseComplement = "CGAT"
    assert umiExtractionFunctions.find_reverse_complement(standardNucleotides_lowerCase) == standardNucleotides_reverseComplement

def test_genomic_function_extract_top_and_bottom_of_sequence():
    seqExtractionLength = 200
    topSeq = "A" * seqExtractionLength
    middleSeq = "T"
    bottomSeq = "C" * seqExtractionLength
    bottomSeq_reverseComplement = umiExtractionFunctions.find_reverse_complement(bottomSeq)
    sequence = topSeq + middleSeq + bottomSeq
    topSeqOutput, bottomSeqOutput = umiExtractionFunctions.extract_top_and_bottom_of_sequence(sequence)
    assert topSeqOutput == topSeq
    assert bottomSeqOutput == bottomSeq_reverseComplement
