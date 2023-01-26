import pytest
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
from general import genomicFunctions
import re

@pytest.fixture
def allIUPAC_nucleotides(): return "ACGTWSRYKMBDHVN"


def test_genomic_function_convert_IUPAC_to_regular_expression(allIUPAC_nucleotides):
    nucleotidesInRegExFormat = "ACGT[AT][CG][AG][CT][GT][AC][CGT][AGT][ACT][ACG][ACGT]"
    assert genomicFunctions.convert_IUPAC_to_regular_expression(allIUPAC_nucleotides) == nucleotidesInRegExFormat

def test_genomic_function_is_only_standard_nucleotide():
    standardNucleotides = "ATCG"
    nonStandardNucleotides = "ATCGR"
    assert genomicFunctions.is_only_standard_nucleotide(standardNucleotides)
    assert not genomicFunctions.is_only_standard_nucleotide(nonStandardNucleotides)

def test_genomic_function_is_IUPAC_nucleotide(allIUPAC_nucleotides):
    assert genomicFunctions.is_IUPAC_nucleotide(allIUPAC_nucleotides)

def test_genomic_function_is_IUPAC_nucleotide_rejects_inappropriate_character(allIUPAC_nucleotides):
    allIUPAC_nucleotides_withError = allIUPAC_nucleotides + "@"
    assert not genomicFunctions.is_IUPAC_nucleotide(allIUPAC_nucleotides_withError)

def test_genomic_function_find_reverse_complement(allIUPAC_nucleotides):
    allIUPAC_nucleotides_reverseComplement = "NBDHVKMRYSWACGT"
    assert genomicFunctions.find_reverse_complement(allIUPAC_nucleotides) == allIUPAC_nucleotides_reverseComplement

def test_genomic_function_find_reverse_complement_non_nucleotide_character_error():
    nonNucleotideCharacter = "@"
    errorOutput = "Provided sequence contains non-nucleotide characters"
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        genomicFunctions.find_reverse_complement(nonNucleotideCharacter)

def test_genomic_function_find_reverse_complement_works_regardless_of_case():
    standardNucleotides_lowerCase = "atcg"
    standardNucleotides_reverseComplement = "CGAT"
    assert genomicFunctions.find_reverse_complement(standardNucleotides_lowerCase) == standardNucleotides_reverseComplement

def test_genomic_function_extract_top_and_bottom_of_sequence():
    seqExtractionLength = 200
    topSeq = "A" * seqExtractionLength
    middleSeq = "T"
    bottomSeq = "C" * seqExtractionLength
    bottomSeq_reverseComplement = genomicFunctions.find_reverse_complement(bottomSeq)
    sequence = topSeq + middleSeq + bottomSeq
    topSeqOutput, bottomSeqOutput = genomicFunctions.extract_top_and_bottom_of_sequence(sequence)
    assert topSeqOutput == topSeq
    assert bottomSeqOutput == bottomSeq_reverseComplement
