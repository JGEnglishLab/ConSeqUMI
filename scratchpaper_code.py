import consensus_maker as cm

seqs = ['AAAAAAAAAAAAAAAAAAAA','TTTTTTTTTTTTTTTTTTTT','CCAAAAAAAAAAAAAAAAAA']

for groups in cm.cluster_longread_consensus_sequences(seqs, threshold = 1/20): print(groups)
