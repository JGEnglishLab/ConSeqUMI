import unittest
import sys, os
sys.path.append(os.path.abspath(os.path.join('..', 'longread_umi_python')))
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import umi_binning as ub

class MyTest(unittest.TestCase):

    def test_cutadapter_adapter_initialization(self):
        forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
        forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
        reverseAdapter1 = 'AATGATACGGCGACCACCGAGATC'
        reverseAdapter2 = 'CGACATCGAGGTGCCAAAC'
        forwardAdapter1_reverseComplement = 'ATCTCGTATGCCGTCTTCTGCTTG'
        forwardAdapter2_reverseComplement = 'CTGAGCCAKRATCRAACYCT'
        reverseAdapter1_reverseComplement = 'GATCTCGGTGGTCGCCGTATCATT'
        reverseAdapter2_reverseComplement = 'GTTTGGCACCTCGATGTCG'

        self.UMIBins = ub.UMIBinner()
        self.UMIBins.set_adapters_for_future_matching(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)
        self.assertEqual(self.UMIBins.forward.front_adapter.sequence, forwardAdapter1)
        self.assertEqual(self.UMIBins.forward.back_adapter.sequence, forwardAdapter2)
        self.assertEqual(self.UMIBins.forward_pair.front_adapter.sequence, reverseAdapter2_reverseComplement)
        self.assertEqual(self.UMIBins.forward_pair.back_adapter.sequence, reverseAdapter1_reverseComplement)
        self.assertEqual(self.UMIBins.reverse.front_adapter.sequence, reverseAdapter1)
        self.assertEqual(self.UMIBins.reverse.back_adapter.sequence, reverseAdapter2)
        self.assertEqual(self.UMIBins.reverse_pair.front_adapter.sequence, forwardAdapter2_reverseComplement)
        self.assertEqual(self.UMIBins.reverse_pair.back_adapter.sequence, forwardAdapter1_reverseComplement)

    def test_adapter_indexer_locator(self):
        self.UMIBins = ub.UMIBinner()
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

        self.UMIBins.set_adapters_for_future_matching(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)

        indices1 = self.UMIBins.get_adapter_indices(exampleForwardSeq)
        self.assertEqual(indices1, (110,100))
        indices2 = self.UMIBins.get_adapter_indices(exampleReverseSeq)
        self.assertEqual(indices2, (110,100))
        indices3 = self.UMIBins.get_adapter_indices(fillerSeq2)
        self.assertEqual(indices3, (-1,-1))

    def test_set_forward_reverse_index(self):
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
        self.UMIBins = ub.UMIBinner()
        self.UMIBins.set_adapters_for_future_matching(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)

        self.UMIBins.set_adapter_indices(indices_array, indel = 5)
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
        self.UMIBins.set_adapter_indices(indices_array2, indel = 5)
        self.assertEqual(self.UMIBins.forward_adapter_start_index, 0)
        self.assertEqual(self.UMIBins.forward_adapter_stop_index, 67)
        self.assertEqual(self.UMIBins.reverse_adapter_start_index, 0)
        self.assertEqual(self.UMIBins.reverse_adapter_stop_index, 66)

    def test_adapter_not_found_assertion(self):

        indices_array = np.array([
            [20, 40],
            [20, 40],
            [20, 40],
            [20, 40],
        ])
        with self.assertRaises(Exception) as context:
            ub.UMIBinner().set_adapter_indices(indices_array, indel = 5)
        self.assertEqual('Adapter matched to fewer than 5 given sequences', str(context.exception))

    def test_sequence_adapter_matching(self):
        #forward read
        fillerSeq1 = 'A'*110
        fillerSeq2 = 'C'*200
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

        exampleErrorSeq = "".join([
                        fillerSeq1,
                        reverseAdapter1, umi2r, 'CGCG', reverseAdapter2,
                        fillerSeq2,
                        forwardAdapter2_reverseComplement, umi1r, forwardAdapter1_reverseComplement +
                        fillerSeq3])

        self.UMIBins = ub.UMIBinner()
        self.UMIBins.set_adapters_for_future_matching(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)
        self.UMIBins.forward_adapter_start_index = 100
        self.UMIBins.forward_adapter_stop_index = 187
        self.UMIBins.reverse_adapter_start_index = 90
        self.UMIBins.reverse_adapter_stop_index = 176
        umiForwardPair = self.UMIBins.extract_umi_from_linked_adapters(exampleForwardSeq)
        umiReversePair = self.UMIBins.extract_umi_from_linked_adapters(exampleReverseSeq)
        umiError = self.UMIBins.extract_umi_from_linked_adapters(exampleErrorSeq)
        self.assertEqual(umiForwardPair, umi1 + umi2)
        self.assertEqual(umiReversePair, umi1 + umi2)
        self.assertIsNone(umiError)

        self.UMIBins.forward_adapter_start_index = 0
        self.UMIBins.forward_adapter_stop_index = 90
        self.UMIBins.reverse_adapter_start_index = 0
        self.UMIBins.reverse_adapter_stop_index = 90

        exampleForwardSeq2 = exampleForwardSeq[110:-100]
        exampleReverseSeq2 = exampleReverseSeq[100:-110]

        umiForwardPair2 = self.UMIBins.extract_umi_from_linked_adapters(exampleForwardSeq2)
        umiReversePair2 = self.UMIBins.extract_umi_from_linked_adapters(exampleReverseSeq2)
        self.assertEqual(umiForwardPair2, umi1 + umi2)
        self.assertEqual(umiReversePair2, umi1 + umi2)

    def test_too_few_umis_identified(self):
        umi_sequences = [
            'AAAAGGCACTAGTGCCCA',
            'TTAAGGCACTAGTGCCCT',
            'AAAAGGCACTAGTGCCCA',
            'TTAAGGCACTAGTGCCCT',
        ]
        with self.assertRaises(Exception) as context:
            ub.UMIBinner().find_consensus_sequences(umi_sequences, 5)
        self.assertEqual('Fewer than 5 UMI sequences match UMI pattern', str(context.exception))

    def test_hamming_distance_matrix(self):
        seqs = ['ATCG','AACG','TCGC', 'TTCC']
        distanceMatrix1 = np.array([0.25, 1, 0.5, 1, 0.25, 0.5]) #0,1; 0,2; 0,3; 1,2; 1,3; 2,3
        distanceMatrix2 = ub.UMIBinner().make_hamming_distance_matrix(seqs)
        self.assertEqual(distanceMatrix2.all(), distanceMatrix1.all())

    def test_consensus_assignment(self):
        forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
        forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
        reverseAdapter1 = 'AATGATACGGCGACCACCGAGATC'
        reverseAdapter2 = 'CGACATCGAGGTGCCAAAC'
        self.UMIBins = ub.UMIBinner()
        self.UMIBins.is_umi_patterned = False
        self.UMIBins.set_adapters_for_future_matching(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)

        consensus_sequences1 = ['TAAAGGCACTAGTGCCCTTGCCGTATAAGGTACTCG','AAGATGCATCGACTCGTCCGCTATTGTGAACGTGCA']
        #consensus_labels1 = [1,2,1,2,1,2,1,2,1,2,1,2,1,2]
        #consensus_labels1 = [0,1,0,1,0,1,0,1,0,1,0,1,0,1]
        consensus_labels1 = [2,1,2,1,2,1,2,1,2,1,2,1,2,1]

        consensus_sequences2, consensus_labels2 = self.UMIBins.identify_consensus_umi_sequences_from_files(['test/tdd_example.fq'], min_clusters=2)

        self.assertSetEqual(set(consensus_sequences2), set(consensus_sequences1))
        self.assertListEqual(list(consensus_labels2), consensus_labels1)

    def test_too_few_clusters_identified(self):
        umi_sequences = [ #4 clusters of 4 each - none should succeed because clusters have < 5
            'AGATCAAGCTTTAGCCGCTGACTGTACGAGTGACTA',
            'CTATCAAGCTTTAGCCGCTGACTGTACGAGTGACTA',
            'CGCTCAAGCTTTAGCCGCTGACTGTACGAGTGACTA',
            'CGAGCAAGCTTTAGCCGCTGACTGTACGAGTGACTA',

            'AATATTAAGGTCAGCTGAACCGCTTAGAGCCGGCTT',
            'CTTATTAAGGTCAGCTGAACCGCTTAGAGCCGGCTT',
            'CACATTAAGGTCAGCTGAACCGCTTAGAGCCGGCTT',
            'CATGTTAAGGTCAGCTGAACCGCTTAGAGCCGGCTT',

            'AGGTGGGTGCGTACAATCTCTCGAACCATATAACGT',
            'CTGTGGGTGCGTACAATCTCTCGAACCATATAACGT',
            'CGCTGGGTGCGTACAATCTCTCGAACCATATAACGT',
            'CGGGGGGTGCGTACAATCTCTCGAACCATATAACGT',

            'ATCGGTAGGGTGCCTACAAGTCCAACTCGTTAGACA',
            'TGCGGTAGGGTGCCTACAAGTCCAACTCGTTAGACA',
            'TTTGGTAGGGTGCCTACAAGTCCAACTCGTTAGACA',
            'TTCCGTAGGGTGCCTACAAGTCCAACTCGTTAGACA',
        ]

        with self.assertRaises(Exception) as context:
            ub.UMIBinner().find_consensus_sequences(umi_sequences, min_clusters=5)
        self.assertEqual('Fewer than 5 clusters contain at least 5 UMI sequences', str(context.exception))


def create_test_file():
    u = ub.UMIBinner()

    fillerSeq1 = 'A'*30
    fillerSeq2 = 'C'*200
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
    ]
    umi_sequences2 = [
        'ATGATGCATCGACTCGTA',
        'AACATGCATCGACTCGTC',
        'AAGGTGCATCGACTCGTA',
        'AAGAAGCATCGACTCGTC',
        'AAGATTCATCGACTCGTC',
        'AAGATGGATCGACTCGTC',
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
    ]
    umi_sequences4 = [
        'AGCACGTTCACAATAGCA',
        'TTCACGTTCACAATAGCG',
        'TGGACGTTCACAATAGCA',
        'TGCCCGTTCACAATAGCG',
        'TGCAAGTTCACAATAGCG',
        'TGCACTTTCACAATAGCG',
        'TGCACGTTCACAATAGCG',
    ]

    sequences = []
    for i in range(7):
        seq1 = ''.join([
            fillerSeq1,
            forwardAdapter1, umi_sequences1[i], forwardAdapter2,
            fillerSeq2,
            rc_reverseAdapter2, u.reverse_complement(umi_sequences3[i]), rc_reverseAdapter1,
            fillerSeq3
        ])
        seq2 = ''.join([
            fillerSeq1,
            forwardAdapter1, umi_sequences2[i], forwardAdapter2,
            fillerSeq2,
            rc_reverseAdapter2, u.reverse_complement(umi_sequences4[i]), rc_reverseAdapter1,
            fillerSeq3
        ])
        if i % 3 == 0:
            seq1 = u.reverse_complement(seq1)
            seq2 = u.reverse_complement(seq2)

        rec1 = SeqRecord(Seq(seq1),id=str(i))
        rec2 = SeqRecord(Seq(seq2),id=str(i+8))
        rec1.letter_annotations["phred_quality"] = [1 for x in seq1]
        rec2.letter_annotations["phred_quality"] = [1 for x in seq2]

        sequences.extend([rec1, rec2])

    errorRec = SeqRecord(Seq('AAAAAAAAAAAAATTTTTTT'), id='error')
    errorRec.letter_annotations["phred_quality"] = [1 for x in 'AAAAAAAAAAAAATTTTTTT']
    sequences.append(errorRec)
    with open("test/tdd_example.fq", "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fastq")

if __name__ == '__main__':
    #create_test_file()
    unittest.main()
