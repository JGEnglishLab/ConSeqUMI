import umi_binning as ub
import pickle
import numpy as np
import os
from collections import defaultdict
'''
forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
reverseAdapter1 = 'AATGATACGGCGACCACCGAGATC'
reverseAdapter2 = 'CGACATCGAGGTGCCAAAC'
#file_name = 'test/data/test_reads.fq'
#file_name = 'test/data/ont_r10_zymo_rrna.fq'
#file_name = 'test/data/pb_ccs_zymo_rrna.fq'
#'''
#'''
forwardAdapter1 = 'GAGTGTGGCTCTTCGGAT'
forwardAdapter2 = 'ATCTCTACGGTGGTCCTAAATAGT'
reverseAdapter1 =  'GTGGGACTGCTGATGACGACTGAT'
reverseAdapter2 = 'GGCGCGTTTTTTTTTTTTTTTTTT'
file_name = 'test/data/fastq_runid_03e0c9c93c94bab9a212849ca3c0e8409f5dc160_0_0.fastq'
#'''


UMIBins = ub.UMIBinner()
UMIBins.set_adapters_for_future_matching(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)

dir = 'test/data/UMI fastq for Caleb/'
files = [dir+file for file in os.listdir(dir)]
#files = [file_name]
print('enter function')
UMIBins.identify_adapter_start_end_indices(files)
umi_sequences = UMIBins.identify_umi_sequences(files)
print(len(umi_sequences))
with open('/Users/calebcranney/Documents/adam_1_seq.txt', 'w') as f:
    for item in umi_sequences:
        f.write("%s\n" % item)
consensus_sequences, consensus_labels = UMIBins.identify_consensus_umi_sequences_from_files(files)
pickle.dump( consensus_sequences, open( "test/data/consensus_sequences.p", "wb" ) )
pickle.dump( consensus_labels, open( "test/data/consensus_labels.p", "wb" ) )
np.set_printoptions(threshold=np.inf)
print(consensus_labels)
conDict = defaultdict(int)
for i in consensus_labels:
    conDict[i] += 1
counts = sorted([conDict[i] for i in conDict], reverse=True)
print(len(consensus_sequences))
print(max(consensus_labels))
print('\n')
for i in counts: print(i)
'''
#consensus = UMIBins.identify_consensus_umi_sequences_from_file(file_name)
#pickle.dump( consensus, open( "test/data/consensus.p", "wb" ) )
consensus = pickle.load(open( "test/data/consensus.p", "rb" ))
umi = pickle.load(open( "test/data/umi_sequences.p", "rb" ))
umi = np.array(umi)
for cluster_id in np.unique(consensus.labels_):
    exemplar = umi[consensus.cluster_centers_indices_[cluster_id]]
    cluster = np.unique(umi[np.nonzero(consensus.labels_==cluster_id)])
    cluster_str = ", \n".join(cluster)
    print(len(cluster))
    print(" - *%s:* \n%s\n\n" % (exemplar, cluster_str))

print(len(consensus.labels_))
print(consensus.labels_)
'''
