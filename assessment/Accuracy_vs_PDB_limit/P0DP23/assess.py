

from Bio import PDB, SeqIO
from Bio.SeqUtils import IUPACData
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as sp
import argparse

parser = argparse.ArgumentParser(description = "Open protein FASTA files")
parser.add_argument('-s', '--suffix',
                        dest = 'suffix', action = 'store',
                        default = None,
                        help = 'suffix of the input file',
                        type = str)

options = parser.parse_args()

suffix = options.suffix

## beta factors from prediction ##
prediction = pd.read_csv("P0DP23_"+suffix+".txt", skiprows=1, sep="\t", header=None)
prediction = prediction[[0, 2]]
prediction.columns = ['position', 'bfactor']
prediction = prediction.assign(key="pred")


## beta factors from pdb ##
parser=PDB.PDBParser()
structure = parser.get_structure("1cll", "1cll.pdb")
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
plt.savefig('comparison_'+suffix+'.png')


plt.figure()
plot = sns.regplot(x="bfactor_x", y="bfactor_y", data = data2)
plot.set_title("Correlation between predicted and observed bfactors")
plot.set_xlabel("observed")
plot.set_ylabel("predicted")
plt.savefig('correlation_'+suffix+'.png')

# Pearson correlation test
x = data2['bfactor_x'].to_numpy()
y = data2['bfactor_y'].to_numpy()
nas = np.logical_or(np.isnan(x), np.isnan(y))
with open("P0DP23_"+suffix+".pearsonr", 'wt') as fd:
    fd.write(suffix+str(sp.pearsonr(x[~nas], y[~nas])))
