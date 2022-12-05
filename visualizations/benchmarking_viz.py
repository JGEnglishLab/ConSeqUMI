import matplotlib.pyplot as plt
import pandas as pd
import sys

benchmarkPath = sys.argv[1]
df = pd.read_csv(benchmarkPath)
df = df[['clusterSize','levenshteinDistance']]
df = df.groupby('clusterSize').mean()
df.reset_index(inplace=True)
print(df)

lines = df.plot.line(x='clusterSize', y='levenshteinDistance')
output = benchmarkPath.split('.')[0] + '_hist.png'
#plt.show()
plt.savefig(output)
