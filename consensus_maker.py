from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import numpy as np
from collections import Counter
import unittest
from Bio import AlignIO
import subprocess

def find_consensus_seq_from_list(seqs, index, umiBinPath):
    records = (SeqRecord(Seq(seq.strip()), str(index)) for index,seq in enumerate(seqs))
    infile = umiBinPath.split('.')[0] + str(index) + 'unaligned.fa'
    outfile = umiBinPath.split('.')[0] + str(index) + 'aligned.fa'
    SeqIO.write(records, infile, 'fasta')

    process = subprocess.Popen(['/usr/local/anaconda3/envs/longread/bin/muscle',
     '-in', infile,
     '-out', outfile])
    stdout, stderr = process.communicate()

    alignment = AlignIO.read(outfile, format='fasta')
    #for i in range(len(seqs)): alignment.add_sequence(str(i), seqs[i])
    info = AlignInfo.SummaryInfo(alignment)
    return info.dumb_consensus(threshold=0.7)

def find_consensus_sequences_from_umi_bins(umiBinPath, seqPath):
    umiBins = pd.read_csv(umiBinPath, sep='\t')
    umiBins.columns = ['umi','count','indices']
    trimmedBins = umiBins[umiBins['count'] > 10]

    with open(seqPath, 'r') as file: sequences = np.array(file.readlines())

    finalSeqs = []
    for index, row in trimmedBins.iterrows():
        indices = row['indices'].split(',')
        indices = [int(i)-1 for i in indices]
        cluster = sequences[indices]
        consensusSequence = find_consensus_seq_from_list(cluster, index, umiBinPath)
        finalSeqs.append([row['umi'],str(consensusSequence).strip('\n')])

    return pd.DataFrame(finalSeqs, columns=['UMIBin','ConsensusSequence'])
