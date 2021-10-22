from Bio.Align import MultipleSeqAlignment, AlignInfo
import pandas as pd
import numpy as np

def find_consensus_seq_from_list(seqs):
    alignment = MultipleSeqAlignment([])
    for i in range(len(seqs)): alignment.add_sequence(str(i), seqs[i])
    info = AlignInfo.SummaryInfo(alignment)
    return info.dumb_consensus()

def find_consensus_sequences_from_umi_bins(umiBinPath, seqPath):
    umiBins = pd.read_csv(umiBinPath, sep='\t')
    umiBins.columns = ['umi','count','indices']
    trimmedBins = umiBins[umiBins['count'] > 4]

    with open(seqPath, 'r') as file: sequences = np.array(file.readlines())

    finalSeqs = []
    for index, row in trimmedBins.iterrows():
        indices = row['indices'].split(',')
        indices = [int(i)-1 for i in indices]
        cluster = sequences[indices]
        consensusSequence = find_consensus_seq_from_list(cluster)
        finalSeqs.append([row['umi'],str(consensusSequence).strip('\n')])
    data = [['TAAAGGCACTAGTGCCCTTGCCGTATAAGGTACTCG','CCTTGTCTTGCGTCGTCCCCTGCGTATGCCTCTCGAGTCGATCATCAGACCTATGCAGGGCGCGGTGTGTAAGGTACCCCCCGTGTCCCGCGTGCCCAGTGCAGAACTCAAGAGCGGAGGGTTCAACAAACCATATAGTAGAACATGCCAGCAGCTTGTATTGGTTCTGAGTATACTTTCGGGTTGAATAGTCGTTGTGA'],
            ['AAGATGCATCGACTCGTCCGCTATTGTGAACGTGCA','TTAGTGCCTTTCCCACGAGCTGTAGACGAAGCTCACTATTCCAGGTCTAAGGCCCAGGGGCGTGAGAATCGTGAACTTGTTAAAGATTTACTCCTTGTTGCCGTCGAAGAATGTTGGTTCTCTTACGCTGATAACAAATGACCATGCTCTTTGAGATCCCAATCCGGCAACAAAACTAAGCGCGGGCCCGCAGGCTCTGC']]

    return pd.DataFrame(finalSeqs, columns=['UMIBin','ConsensusSequence'])
