

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
prediction = prediction.assign(key="pred")


## beta factors from pdb ##
parser=PDB.PDBParser()
structure = parser.get_structure("1x1r", "1x1r.pdb")
structure = structure[0]["A"]

### NORMALIZE BFACTORS
bfactors = list()
for residue in structure.get_residues():
    if "CA" in residue:
        bfactors.append(residue['CA'].get_bfactor())
bf_mean = np.mean(bfactors)
bf_sd = np.std(bfactors)
for residue in structure.get_residues():
    if "CA" in residue:
        residue['CA'].set_bfactor((residue['CA'].get_bfactor() - bf_mean)/bf_sd)

### get bfactors
obs_data = []
for residue in structure.get_residues():
    if "CA" in residue:
        pos = residue.get_id()[1]
        b = residue['CA'].get_bfactor()
        obs_data.append([pos, b])

observed = pd.DataFrame(obs_data, columns=['position', 'bfactor'])
observed = observed.assign(key="obs")

### plot data ###
data = pd.concat([prediction, observed], axis=0, ignore_index=True)

data2 = observed.merge(prediction, on="position")

plt.figure()
plot = sns.lineplot(x='position', y='bfactor', hue='key', data=data)
plot.set_title("Comparison between predicted and pdb")
plot.set_xlabel("Residue number")
plot.set_ylabel("B-factor")
plt.savefig('comparison.png')


plt.figure()
plot = sns.regplot(x="bfactor_x", y="bfactor_y", data = data2)
plot.set_title("Correlation between predicted and observed bfactors")
plot.set_xlabel("observed")
plot.set_ylabel("predicted")
plt.savefig('correlation.png')

# Pearson correlation test
x = data2['bfactor_x'].to_numpy()
y = data2['bfactor_y'].to_numpy()
nas = np.logical_or(np.isnan(x), np.isnan(y))
print(sp.pearsonr(x[~nas], y[~nas]))
