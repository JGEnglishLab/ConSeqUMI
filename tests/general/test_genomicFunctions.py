import pytest
import sys
sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
from general import genomicFunctions

def test_genomic_function_convert_IUPAC_to_regular_expression():
    allIUPAC_nucleotides = "ACGTWSRYKMBDHVN"
    nucleotidesInRegExFormat = "ACGT[AT][CG][AG][CT][GT][AC][CGT][AGT][ACT][ACG][ACGT]"
    assert genomicFunctions.convert_IUPAC_to_regular_expression(allIUPAC_nucleotides) == nucleotidesInRegExFormat

def test_genomic_function_find_reverse_complement():
    allIUPAC_nucleotides = "ACGTWSRYKMBDHVN"
    allIUPAC_nucleotides_reverseComplement = "NBDHVKMRYSWACGT"

    assert genomicFunctions.find_reverse_complement(allIUPAC_nucleotides) == allIUPAC_nucleotides_reverseComplement

def test_genomic_function_find_reverse_complement_accounts_for_only_standard_nucleotides_when_specified():
    standardNucleotides = "ATCG"
    standardNucleotides_reverseComplement = "CGAT"
    nonStandardNucleotides = "ATCGR"
    nonStandardNucleotides_reverseComplement = "YCGAT"
    assert genomicFunctions.find_reverse_complement(standardNucleotides, isOnlyStandardNucleotide=True) == standardNucleotides_reverseComplement
    assert genomicFunctions.find_reverse_complement(nonStandardNucleotides, isOnlyStandardNucleotide=False) == nonStandardNucleotides_reverseComplement
    with pytest.raises(TypeError):
        genomicFunctions.find_reverse_complement(nonStandardNucleotides, isOnlyStandardNucleotide=True)
