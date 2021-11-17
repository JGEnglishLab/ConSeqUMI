from Bio import SeqIO
import consensus_maker as cm

consensusFile = 'test/data/output/consensus.fasta'
seqs = [str(record.seq) for record in SeqIO.parse(consensusFile, "fasta")]

labels = [x for x in cm.cluster_longread_consensus_sequences(seqs, dendrogramFile='test/data/output/dendrogram.png')]
