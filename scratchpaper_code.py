import umi_binning as ub

forwardAdapter1 = 'CAAGCAGAAGACGGCATACGAGAT'
forwardAdapter2 = 'AGRGTTYGATYMTGGCTCAG'
reverseAdapter1 = 'AATGATACGGCGACCACCGAGATC'
reverseAdapter2 = 'CGACATCGAGGTGCCAAAC'
UMIBins = ub.UMIBinner()
UMIBins.set_adapters_for_future_matching(forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2)

consensus2 = UMIBins.identify_consensus_umi_sequences_from_file('test/test_reads.fq')
print(consensus2)
