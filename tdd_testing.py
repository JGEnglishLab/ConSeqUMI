import unittest
import umi_binning as ub

class MyTest(unittest.TestCase):

    def test_reverse_adapter_code(self): self.assertEqual(ub.make_reverse_complementary_adapter_sequence('CAAGCAGAAGACGGCATACGAGAT'), 'ATCTCGTATGCCGTCTTCTGCTTG')
    


#don't alter - this allows running 'python TDD_template.py' to run the test cases
if __name__ == '__main__':
    unittest.main()
