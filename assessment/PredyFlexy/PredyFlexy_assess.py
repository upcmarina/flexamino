

from Bio import PDB, SeqIO
from Bio.SeqUtils import IUPACData
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as sp



## beta factors from prediction ##
prediction = pd.read_csv("O08989_pred_out.txt", skiprows=1, sep="\t", header=None)
prediction = prediction[[0, 2]]
prediction.columns = ['position', 'bfactor']
prediction = prediction.assign(key="FlexAmino")


## beta factors from PredyFlexy ##
PredyFlexy = pd.read_csv("O08989.predyflexy", skiprows=1, sep="\s+",header=None)
observed = PredyFlexy[[0, 9]]
observed.columns = ['position', 'bfactor']
observed = observed.assign(key="PredyFlexy")

observed.loc[observed['bfactor'] == -9.999, 'bfactor'] = np.nan

### plot data ###
data = pd.concat([prediction, observed], axis=0, ignore_index=True)

data2 = observed.merge(prediction, on="position")

plt.figure()
plot = sns.lineplot(x='position', y='bfactor', hue='key', data=data)
plot.set_title("Comparison between FlexAmino and PredyFlexy scores")
plot.set_xlabel("Residue number")
plot.set_ylabel("B-factor")
plt.savefig('comparison.png')


plt.figure()
plot = sns.regplot(x="bfactor_x", y="bfactor_y", data = data2)
plot.set_title("Correlation between FlexAmino and PredyFlexy bfactors")
plot.set_xlabel("PredyFlexy")
plot.set_ylabel("FlexAmino")
plt.savefig('correlation.png')

### Correlation test
x = observed['bfactor'].to_numpy()
y = prediction['bfactor'].to_numpy()
nas = np.logical_or(np.isnan(x), np.isnan(y))
print(sp.pearsonr(x[~nas], y[~nas]))
