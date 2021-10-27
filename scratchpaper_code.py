import matplotlib.pyplot as plt
import pandas as pd

umiBinPath = '~/Documents/test_4_output_s.txt'
umiBins = pd.read_csv(umiBinPath, sep='\t')
umiBins = umiBins[umiBins.iloc[:, 1] > 4]

plt.hist(umiBins.iloc[:, 1])
plt.ylabel('Count')
plt.xlabel('Cluster Size')
plt.savefig('test_4_output_s_gt4.png')
