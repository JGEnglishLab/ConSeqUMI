import unittest
import sys, os
sys.path.append(os.path.abspath(os.path.join('..', 'longread_umi_python')))
from Bio import SeqIO
from Bio.Seq import Seq
import difflib
from collections import defaultdict

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

        ub.UMIBins = ub.UMIBinner()
        ub.UMIBins.set_adapters_for_future_matching(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)
        self.assertEqual(ub.UMIBins.forward.front_adapter.sequence, forwardAdapter1)
        self.assertEqual(ub.UMIBins.forward.back_adapter.sequence, forwardAdapter2)
        self.assertEqual(ub.UMIBins.forward_pair.front_adapter.sequence, reverseAdapter2_reverseComplement)
        self.assertEqual(ub.UMIBins.forward_pair.back_adapter.sequence, reverseAdapter1_reverseComplement)
        self.assertEqual(ub.UMIBins.reverse.front_adapter.sequence, reverseAdapter1)
        self.assertEqual(ub.UMIBins.reverse.back_adapter.sequence, reverseAdapter2)
        self.assertEqual(ub.UMIBins.reverse_pair.front_adapter.sequence, forwardAdapter2_reverseComplement)
        self.assertEqual(ub.UMIBins.reverse_pair.back_adapter.sequence, forwardAdapter1_reverseComplement)

    def test_sequence_adapter_matching(self):
        #forward read
        fillerSeq1 = 'A'*120
        fillerSeq2 = 'C'*200
        fillerSeq3 = 'T'*120
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
                        fillerSeq1,
                        reverseAdapter1, umi2r, reverseAdapter2,
                        fillerSeq2,
                        forwardAdapter2_reverseComplement, umi1r, forwardAdapter1_reverseComplement +
                        fillerSeq3])

        ub.UMIBins = ub.UMIBinner()
        ub.UMIBins.set_adapters_for_future_matching(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)
        umiForwardPair = ub.UMIBins.match_sequence_to_adapters(exampleForwardSeq)
        umiReversePair = ub.UMIBins.match_sequence_to_adapters(exampleReverseSeq)
        self.assertEqual(umiForwardPair, umi1 + umi2)
        self.assertEqual(umiReversePair, umi1 + umi2)
        #self.assertTrue(ub.UMIBins.umi_pattern.match(umiForwardPair))

    def test_sequence_adapter_matching(self):
        forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
        forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
        reverseAdapter1 = 'AATGATACGGCGACCACCGAGATC'
        reverseAdapter2 = 'CGACATCGAGGTGCCAAAC'
        testFilePath = 'test/test_reads.fq'
        ansFilePath = 'test/umi12f.fa'
        ub.UMIBins = ub.UMIBinner()
        ub.UMIBins.set_adapters_for_future_matching(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)
        resultGenerator = ub.UMIBins.extract_umi_pairs_from_file(testFilePath)

        ansDict = defaultdict(int)
        resultDict = defaultdict(int)
        results = [result for result in resultGenerator if result != None]
        #for result in resultGenerator: print('res: '+result)#resultDict[result] += 1
        with open(ansFilePath) as ansFile:
            ansFileParser = SeqIO.parse(ansFile,"fasta")
            #for ans in ansFileParser: print('ans: '+str(ans.seq))#ansDict[str(ans.seq)] += 1
            answers = [str(ans.seq) for ans in ansFileParser]

        print(len(results))
        print(len(answers))
        for i in range(len(results)):
            print(results[i])
            print(answers[i])
            print('\n')





    #test umi are length=18

if __name__ == '__main__':
    unittest.main()
