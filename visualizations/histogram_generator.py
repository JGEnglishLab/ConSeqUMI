import matplotlib.pyplot as plt
import pandas as pd
import sys

umiBinPath = sys.argv[1]
umiBins = pd.read_csv(umiBinPath, sep="\t")
# umiBins = umiBins[umiBins.iloc[:, 1] > 4]


plt.hist(umiBins.iloc[:, 1])
plt.ylabel("Count")
plt.xlabel("Cluster Size")
output = umiBinPath.split(".")[0] + "_hist.png"
plt.savefig(output)
