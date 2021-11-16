import unittest
import sys, os
sys.path.append(os.path.abspath(os.path.join('..', 'longread_umi_python')))
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import umi_extractor as ue
import consensus_maker as cm
import random
import filecmp
import pandas as pd

class MyTest(unittest.TestCase):

    def test_ue_cutadapter_adapter_initialization(self):
        forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
        forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
        reverseAdapter1 = 'AATGATACGGCGACCACCGAGATC'
        reverseAdapter2 = 'CGACATCGAGGTGCCAAAC'
        forwardAdapter1_reverseComplement = 'ATCTCGTATGCCGTCTTCTGCTTG'
        forwardAdapter2_reverseComplement = 'CTGAGCCAKRATCRAACYCT'
        reverseAdapter1_reverseComplement = 'GATCTCGGTGGTCGCCGTATCATT'
        reverseAdapter2_reverseComplement = 'GTTTGGCACCTCGATGTCG'

        self.UMIBins = ue.UMIExtractor()
        self.UMIBins.set_universal_front_and_reverse_linked_adapters(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)
        self.assertEqual(self.UMIBins.forward.front_adapter.sequence, forwardAdapter1)
        self.assertEqual(self.UMIBins.forward.back_adapter.sequence, forwardAdapter2)
        self.assertEqual(self.UMIBins.forward_pair.front_adapter.sequence, reverseAdapter2_reverseComplement)
        self.assertEqual(self.UMIBins.forward_pair.back_adapter.sequence, reverseAdapter1_reverseComplement)
        self.assertEqual(self.UMIBins.reverse.front_adapter.sequence, reverseAdapter1)
        self.assertEqual(self.UMIBins.reverse.back_adapter.sequence, reverseAdapter2)
        self.assertEqual(self.UMIBins.reverse_pair.front_adapter.sequence, forwardAdapter2_reverseComplement)
        self.assertEqual(self.UMIBins.reverse_pair.back_adapter.sequence, forwardAdapter1_reverseComplement)

    def test_ue_adapter_index_locator(self):
        self.UMIBins = ue.UMIExtractor()
        fillerSeq1 = 'A'*110
        fillerSeq2 = 'C'*400
        fillerSeq3 = 'T'*100
        forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
        forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
        reverseAdapter1 = 'AATGATACGGCGACCACCGAGATC'
        reverseAdapter2 = 'CGACATCGAGGTGCCAAAC'
        forwardAdapter1_reverseComplement = 'ATCTCGTATGCCGTCTTCTGCTTG'
        forwardAdapter2_reverseComplement = 'CTGAGCCAKRATCRAACYCT'
        reverseAdapter1_reverseComplement = 'GATCTCGGTGGTCGCCGTATCATT'
        reverseAdapter2_reverseComplement = 'GTTTGGCACCTCGATGTCG'
        umi1 = 'ATATGTACTATTGTATAC'
        umi2 = 'CAATGTGTCGGCATAGGG'
        umi1r = 'GTATACAATAGTACATAT'
        umi2r = 'CCCTATGCCGACACATTG'

        exampleForwardSeq = "".join([
                        fillerSeq1,
                        forwardAdapter1, umi1, forwardAdapter2,
                        fillerSeq2,
                        reverseAdapter2_reverseComplement, umi2, reverseAdapter1_reverseComplement +
                        fillerSeq3])

        exampleReverseSeq = "".join([
                        fillerSeq3,
                        reverseAdapter1, umi2r, reverseAdapter2,
                        fillerSeq2,
                        forwardAdapter2_reverseComplement, umi1r, forwardAdapter1_reverseComplement +
                        fillerSeq1])

        exampleForwardSeq2 = exampleForwardSeq[110:]
        exampleReverseSeq2 = exampleForwardSeq[:-100]

        self.UMIBins.set_universal_front_and_reverse_linked_adapters(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)

        indices1 = self.UMIBins.get_front_and_reverse_adapter_start_indices(exampleForwardSeq)
        self.assertEqual(indices1, (110,100))
        indices2 = self.UMIBins.get_front_and_reverse_adapter_start_indices(exampleReverseSeq)
        self.assertEqual(indices2, (110,100))
        indices3 = self.UMIBins.get_front_and_reverse_adapter_start_indices(fillerSeq2)
        self.assertEqual(indices3, (-1,-1))

    def test_ue_set_universal_forward_reverse_adapter_indices(self):
        forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
        forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
        reverseAdapter1 = 'AATGATACGGCGACCACCGAGATC'
        reverseAdapter2 = 'CGACATCGAGGTGCCAAAC'

        indices_array = np.array([
            [20, 40],
            [20, 40],
            [20, 40],
            [20, 40],
            [20, 40],
            [20, 40],
            [20, 40],
            [20, 40],
            [20, 40],
            [20, 40],
            [20, 40],
            [20, 16],
            [20, 44],
            [20, 2000],
            [18, 40],
            [22, 40],
            [1000, 40],
        ])
        self.UMIBins = ue.UMIExtractor()
        self.UMIBins.set_universal_front_and_reverse_linked_adapters(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)

        self.UMIBins.set_universal_forward_reverse_adapter_indices(indices_array, indel = 5)
        self.assertEqual(self.UMIBins.forward_adapter_start_index, 15)
        self.assertEqual(self.UMIBins.forward_adapter_stop_index, 87)
        self.assertEqual(self.UMIBins.reverse_adapter_start_index, 35)
        self.assertEqual(self.UMIBins.reverse_adapter_stop_index, 106)

        indices_array2 = np.array([
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 16],
            [0, 44],
            [0, 2000],
            [18, 0],
            [22, 0],
            [1000, 0],
        ])
        self.UMIBins.set_universal_forward_reverse_adapter_indices(indices_array2, indel = 5)
        self.assertEqual(self.UMIBins.forward_adapter_start_index, 0)
        self.assertEqual(self.UMIBins.forward_adapter_stop_index, 67)
        self.assertEqual(self.UMIBins.reverse_adapter_start_index, 0)
        self.assertEqual(self.UMIBins.reverse_adapter_stop_index, 66)

    def test_ue_adapter_not_found_assertion(self):

        indices_array = np.array([
            [20, 40],
            [20, 40],
            [20, 40],
            [20, 40],
        ])
        with self.assertRaises(Exception) as context:
            ue.UMIExtractor().set_universal_forward_reverse_adapter_indices(indices_array, indel = 5)
        self.assertEqual('Adapter matched to fewer than 5 given sequences', str(context.exception))

    def test_ue_umi_sequence_extracted_from_adapters_and_matched(self):
        #forward read
        fillerSeq1 = 'A'*110
        random.seed(0)
        consensusSequence = ''.join(random.choices('ATGC', k=200))
        fillerSeq3 = 'T'*100
        forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
        forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
        reverseAdapter1 = 'AATGATACGGCGACCACCGAGATC'
        reverseAdapter2 = 'CGACATCGAGGTGCCAAAC'
        forwardAdapter1_reverseComplement = 'ATCTCGTATGCCGTCTTCTGCTTG'
        forwardAdapter2_reverseComplement = 'CTGAGCCAKRATCRAACYCT'
        reverseAdapter1_reverseComplement = 'GATCTCGGTGGTCGCCGTATCATT'
        reverseAdapter2_reverseComplement = 'GTTTGGCACCTCGATGTCG'
        umi1 = 'ATATGTACTATTGTATAC'
        umi2 = 'CAATGTGTCGGCATAGGG'
        umi1r = 'GTATACAATAGTACATAT'
        umi2r = 'CCCTATGCCGACACATTG'

        exampleForwardSeq = "".join([
                        fillerSeq1,
                        forwardAdapter1, umi1, forwardAdapter2,
                        consensusSequence,
                        reverseAdapter2_reverseComplement, umi2, reverseAdapter1_reverseComplement +
                        fillerSeq3])

        exampleReverseSeq = "".join([
                        fillerSeq3,
                        reverseAdapter1, umi2r, reverseAdapter2,
                        ue.UMIExtractor().reverse_complement(consensusSequence),
                        forwardAdapter2_reverseComplement, umi1r, forwardAdapter1_reverseComplement +
                        fillerSeq1])

        exampleErrorSeq = "".join([
                        fillerSeq1,
                        reverseAdapter1, umi2r, 'CGC', reverseAdapter2,
                        consensusSequence,
                        forwardAdapter2_reverseComplement, umi1r, forwardAdapter1_reverseComplement +
                        fillerSeq3])

        exampleForwardSeq = SeqRecord(Seq(exampleForwardSeq),id='1')
        exampleReverseSeq = SeqRecord(Seq(exampleReverseSeq),id='2')
        exampleErrorSeq = SeqRecord(Seq(exampleErrorSeq),id='3')

        self.UMIBins = ue.UMIExtractor()
        self.UMIBins.set_universal_front_and_reverse_linked_adapters(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)
        self.UMIBins.forward_adapter_start_index = 100
        self.UMIBins.forward_adapter_stop_index = 187
        self.UMIBins.reverse_adapter_start_index = 90
        self.UMIBins.reverse_adapter_stop_index = 176
        umiForwardPair = self.UMIBins.extract_umi_and_sequence_from_linked_adapters(exampleForwardSeq)
        umiReversePair = self.UMIBins.extract_umi_and_sequence_from_linked_adapters(exampleReverseSeq)
        umiError = self.UMIBins.extract_umi_and_sequence_from_linked_adapters(exampleErrorSeq)
        self.assertEqual(umiForwardPair[0], umi1)
        self.assertEqual(umiForwardPair[1], umi2)
        self.assertEqual(umiForwardPair[2].seq, consensusSequence)
        self.assertEqual(umiReversePair[0], umi1)
        self.assertEqual(umiReversePair[1], umi2)
        self.assertEqual(umiReversePair[2].seq, consensusSequence)
        self.assertIsNone(umiError)

        self.UMIBins.forward_adapter_start_index = 0
        self.UMIBins.forward_adapter_stop_index = 90
        self.UMIBins.reverse_adapter_start_index = 0
        self.UMIBins.reverse_adapter_stop_index = 90

        exampleForwardSeq2 = exampleForwardSeq[110:-100]
        exampleReverseSeq2 = exampleReverseSeq[100:-110]

        umiForwardPairAtZeroIndex = self.UMIBins.extract_umi_and_sequence_from_linked_adapters(exampleForwardSeq2)
        umiReversePairAtZeroIndex = self.UMIBins.extract_umi_and_sequence_from_linked_adapters(exampleReverseSeq2)
        self.assertEqual(umiForwardPairAtZeroIndex[0], umi1)
        self.assertEqual(umiForwardPairAtZeroIndex[1], umi2)
        self.assertEqual(umiForwardPairAtZeroIndex[2].seq, consensusSequence)
        self.assertEqual(umiReversePairAtZeroIndex[0], umi1)
        self.assertEqual(umiReversePairAtZeroIndex[1], umi2)
        self.assertEqual(umiReversePairAtZeroIndex[2].seq, consensusSequence)

    def test_ue_umi_and_sequence_extraction_from_files(self):
        forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
        forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
        reverseAdapter1 = 'AATGATACGGCGACCACCGAGATC'
        reverseAdapter2 = 'CGACATCGAGGTGCCAAAC'
        self.UMIBins = ue.UMIExtractor()
        self.UMIBins.set_universal_front_and_reverse_linked_adapters(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)


        files = ['test/tdd/tdd_example.fq']
        self.UMIBins.identify_and_set_front_and_reverse_adapter_start_indices_from_file(files)
        error = self.UMIBins.extract_umi_and_sequences_from_files(files, 'test/tdd/output')


        self.assertTrue(filecmp.cmp('test/tdd/output/umi1.txt', 'test/tdd/output/tdd_example_umi1_expected.txt', shallow=False))
        filecmp.clear_cache()
        self.assertTrue(filecmp.cmp('test/tdd/output/umi2.txt', 'test/tdd/output/tdd_example_umi2_expected.txt', shallow=False))
        filecmp.clear_cache()
        self.assertTrue(filecmp.cmp('test/tdd/output/seq.fq', 'test/tdd/output/tdd_example_seq_expected.fq', shallow=False))
        filecmp.clear_cache()
        self.assertEqual(error, 1)

    def test_cm_remove_chimeras_from_umi_pairs(self):
        starcode1Path = 'test/tdd/output/tdd_example_starcode1_default_expected_w_chimera.txt'
        starcode2Path = 'test/tdd/output/tdd_example_starcode2_default_expected_w_chimera.txt'
        output = 'test/tdd/output/starcode_without_chimera.txt'

        cm.remove_chimeras_from_umi_pairs(starcode1Path, starcode2Path, output)
        self.assertTrue(filecmp.cmp('test/tdd/output/starcode_without_chimera.txt', 'test/tdd/output/tdd_example_starcode_default_expected_without_chimera.txt', shallow=False))
        filecmp.clear_cache()

    def test_cm_bin_sequences_by_umi_pair(self):
        cm.bin_sequences_by_umi_pair('test/tdd/output/seq.fq','test/tdd/output/starcode_without_chimera.txt')

        self.assertTrue(filecmp.cmp('test/tdd/output/seq_bin0.fq', 'test/tdd/output/tdd_example_seq_expected_bin0.fq', shallow=False))
        filecmp.clear_cache()
        self.assertTrue(filecmp.cmp('test/tdd/output/seq_bin1.fq', 'test/tdd/output/tdd_example_seq_expected_bin1.fq', shallow=False))
        filecmp.clear_cache()

    def test_cm_hamming_distance_matrix(self):
        seqs = ['ATCG','AACG','TCGC', 'TTCC']
        distanceMatrix1 = np.array([0.25, 1, 0.5, 1, 0.25, 0.5]) #0,1; 0,2; 0,3; 1,2; 1,3; 2,3
        distanceMatrix2 = cm.make_hamming_distance_matrix(seqs)
        self.assertEqual(distanceMatrix2.all(), distanceMatrix1.all())

    def test_cm_cluster_longread_consensus_sequences(self):
        seqs = ['AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
                'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
                'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATA',
                'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTTTT']
        expectedResults = np.array([1,2,1,2])
        expectedGroups = np.array([[0,2],[1,3]])
        groups = np.array([group for group in cm.cluster_longread_consensus_sequences(seqs)])
        self.assertEqual(groups.all(),expectedGroups.all())



def strip_non_standard_nucleotide_values(str1, str2):
    if len(str1) != len(str2): raise Exception('consensus sequences have different lengths')
    for i in range(len(str1)-1,-1,-1):
        if str1[i] not in ['A', 'T', 'G', 'C']:
            str1 = str1[:i] + str1[i+1:]
            str2 = str2[:i] + str2[i+1:]
    return str1, str2

def create_imperfect_seq_list(seq, numSeqs, numErrors):
    seqs = [seq for i in range(numSeqs)]
    for i in range(len(seqs)):
        diffIndices = list(set([random.randint(0,len(seq)-1) for j in range(numErrors)]))
        diffNucs = ''.join(random.choices('ATGC', k=len(diffIndices)))
        for j in range(len(diffIndices)):
            k = diffIndices[j]
            seqs[i] = seqs[i][:k]+diffNucs[j]+seqs[i][k+1:]
    return seqs

def create_test_fastq_and_txt_files():
    u = ue.UMIExtractor()
    fillerSeq1 = 'A'*30
    random.seed(0)
    origConsensusSequence1 = ''.join(random.choices('ATGC', k=200))
    origConsensusSequence2 = ''.join(random.choices('ATGC', k=200))
    fillerSeq3 = 'T'*110
    forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
    forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
    rc_reverseAdapter1 = u.reverse_complement('AATGATACGGCGACCACCGAGATC')
    rc_reverseAdapter2 = u.reverse_complement('CGACATCGAGGTGCCAAAC')

    consensus_sequences = ['TAAAGGCACTAGTGCCCTTGCCGTATAAGGTACTCG','AAGATGCATCGACTCGTCCGCTATTGTGAACGTGCA']
    # num errors per consensus sequence: [2, 1, 2, 1, 1, 1, 0]
    umi_sequences1 = [
        'AAAAGGCACTAGTGCCCA',
        'TTAAGGCACTAGTGCCCT',
        'TACAGGCACTAGTGCCCA',
        'TAAGGGCACTAGTGCCCT',
        'TAAAAGCACTAGTGCCCT',
        'TAAAGTCACTAGTGCCCT',
        'TAAAGGCACTAGTGCCCT',
        'TAAAGGCACTAGTGCCCT',
        'TAAAGGCACTAGTGCCCT',
        'TAAAGGCACTAGTGCCCT',
        'TAAAGGCACTAGTGCCCT',
        'TAAAGGCACTAGTGCCCT',
        'TAAAGGCACTAGTGCCCT',
        'TAAAGGCACTAGTGCCCT',
    ]
    umi_sequences2 = [
        'ATGATGCATCGACTCGTA',
        'AACATGCATCGACTCGTC',
        'AAGGTGCATCGACTCGTA',
        'AAGAAGCATCGACTCGTC',
        'AAGATTCATCGACTCGTC',
        'AAGATGGATCGACTCGTC',
        'AAGATGCATCGACTCGTC',
        'AAGATGCATCGACTCGTC',
        'AAGATGCATCGACTCGTC',
        'AAGATGCATCGACTCGTC',
        'AAGATGCATCGACTCGTC',
        'AAGATGCATCGACTCGTC',
        'AAGATGCATCGACTCGTC',
        'AAGATGCATCGACTCGTC',
    ]
    umi_sequences3 = [
        'AGAGTACCTTATACGGCT',
        'CTAGTACCTTATACGGCA',
        'CGGGTACCTTATACGGCT',
        'CGACTACCTTATACGGCA',
        'CGAGAACCTTATACGGCA',
        'CGAGTTCCTTATACGGCA',
        'CGAGTACCTTATACGGCA',
        'CGAGTACCTTATACGGCA',
        'CGAGTACCTTATACGGCA',
        'CGAGTACCTTATACGGCA',
        'CGAGTACCTTATACGGCA',
        'CGAGTACCTTATACGGCA',
        'CGAGTACCTTATACGGCA',
        'CGAGTACCTTATACGGCA',
    ]
    umi_sequences4 = [
        'AGCACGTTCACAATAGCA',
        'TTCACGTTCACAATAGCG',
        'TGGACGTTCACAATAGCA',
        'TGCCCGTTCACAATAGCG',
        'TGCAAGTTCACAATAGCG',
        'TGCACTTTCACAATAGCG',
        'TGCACGTTCACAATAGCG',
        'TGCACGTTCACAATAGCG',
        'TGCACGTTCACAATAGCG',
        'TGCACGTTCACAATAGCG',
        'TGCACGTTCACAATAGCG',
        'TGCACGTTCACAATAGCG',
        'TGCACGTTCACAATAGCG',
        'TGCACGTTCACAATAGCG',
    ]

    consensusSequence1List = create_imperfect_seq_list(origConsensusSequence1, 14, 20)
    consensusSequence2List = create_imperfect_seq_list(origConsensusSequence2, 14, 20)

    sequences = []
    for i in range(len(umi_sequences1)):
        seq1 = ''.join([
            fillerSeq1,
            forwardAdapter1, umi_sequences1[i], forwardAdapter2,
            consensusSequence1List[i],
            rc_reverseAdapter2, u.reverse_complement(umi_sequences3[i]), rc_reverseAdapter1,
            fillerSeq3
        ])
        seq2 = ''.join([
            fillerSeq1,
            forwardAdapter1, umi_sequences2[i], forwardAdapter2,
            consensusSequence2List[i],
            rc_reverseAdapter2, u.reverse_complement(umi_sequences4[i]), rc_reverseAdapter1,
            fillerSeq3
        ])
        if i % 3 == 0:
            seq1 = u.reverse_complement(seq1)
            seq2 = u.reverse_complement(seq2)

        rec1 = SeqRecord(Seq(seq1),id=str(i) + origConsensusSequence1, description='testConsensusSequence1')
        rec2 = SeqRecord(Seq(seq2),id=str(i) + origConsensusSequence2, description='testConsensusSequence2')
        rec1.letter_annotations["phred_quality"] = [1 for x in seq1]
        rec2.letter_annotations["phred_quality"] = [1 for x in seq2]

        sequences.extend([rec1, rec2])

    errorRec = SeqRecord(Seq('AAAAAAAAAAAAATTTTTTT'), id='error')
    errorRec.letter_annotations["phred_quality"] = [1 for x in 'AAAAAAAAAAAAATTTTTTT']
    sequences.append(errorRec)
    with open("test/tdd/tdd_example.fq", "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fastq")

    expectedOutputUMI1 = []
    expectedOutputUMI2 = []
    expectedOutputSeq1 = []
    expectedOutputSeq = []
    expectedOutputError = []
    for i in range(len(umi_sequences1)):
        expectedOutputUMI1.append(umi_sequences1[i])
        expectedOutputUMI2.append(u.reverse_complement(umi_sequences3[i]))
        expectedOutputUMI1.append(umi_sequences2[i])
        expectedOutputUMI2.append(u.reverse_complement(umi_sequences4[i]))

        expectedOutputSeq1.append(consensusSequence1List[i])
        expectedOutputSeq1.append(consensusSequence2List[i])

        seqRec1 = SeqRecord(Seq(consensusSequence1List[i]),id=str(i) + origConsensusSequence1, description='testConsensusSequence1')
        seqRec2 = SeqRecord(Seq(consensusSequence2List[i]),id=str(i) + origConsensusSequence2, description='testConsensusSequence2')
        seqRec1.letter_annotations["phred_quality"] = [1 for x in consensusSequence1List[i]]
        seqRec2.letter_annotations["phred_quality"] = [1 for x in consensusSequence2List[i]]

        expectedOutputSeq.append(seqRec1)
        expectedOutputSeq.append(seqRec2)

    expectedOutputError.append('AAAAAAAAAAAAATTTTTTT')
    with open("test/tdd/output/tdd_example_umi1_expected.txt", "w") as output:
        for seq in expectedOutputUMI1: output.write(seq + '\n')
    with open("test/tdd/output/tdd_example_umi2_expected.txt", "w") as output:
        for seq in expectedOutputUMI2: output.write(seq + '\n')
    with open("test/tdd/output/tdd_example_seq_expected.txt", "w") as output:
        for seq in expectedOutputSeq1: output.write(seq + '\n')
    with open("test/tdd/output/tdd_example_error_expected.txt", "w") as output:
        for seq in expectedOutputError: output.write(seq)
    with open("test/tdd/output/tdd_example_seq_expected.fq", "w") as output_handle:
        SeqIO.write(expectedOutputSeq, output_handle, "fastq")

if __name__ == '__main__':
    #create_test_fastq_and_txt_files()
    unittest.main()
