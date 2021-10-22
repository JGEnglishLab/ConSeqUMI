import umi_extractor as ue
import pickle
import numpy as np
import os
from collections import defaultdict
#'''
forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
reverseAdapter1 = 'AATGATACGGCGACCACCGAGATC'
reverseAdapter2 = 'CGACATCGAGGTGCCAAAC'
file_name = 'test/data/test_reads.fq'
#file_name = 'test/data/ont_r10_zymo_rrna.fq'
#file_name = 'test/data/pb_ccs_zymo_rrna.fq'
#'''
'''
forwardAdapter1 = 'GAGTGTGGCTCTTCGGAT'
forwardAdapter2 = 'ATCTCTACGGTGGTCCTAAATAGT'
reverseAdapter1 =  'GTGGGACTGCTGATGACGACTGAT'
reverseAdapter2 = 'GGCGCGTTTTTTTTTTTTTTTTTT'
file_name = 'test/data/fastq_runid_03e0c9c93c94bab9a212849ca3c0e8409f5dc160_0_0.fastq'
#'''



#dir = 'test/data/UMI fastq for Caleb/'
#files = [dir+file for file in os.listdir(dir)]
files = [file_name]
outputHeader = file_name.split('.')[0]

UMIBins = ue.UMIExtractor()
UMIBins.set_universal_front_and_reverse_linked_adapters(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)
UMIBins.identify_and_set_front_and_reverse_adapter_start_indices_from_file(files)
excludedCount = UMIBins.extract_umi_and_sequences_from_files(files, outputHeader)
print(excludedCount)
