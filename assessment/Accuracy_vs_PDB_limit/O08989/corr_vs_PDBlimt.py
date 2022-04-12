import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


data = pd.read_csv("pearsonR_vs_pbd_lim.txt", skiprows=0, sep=" ", header=None)
data.columns = ['PDB_limit', 'PearsonR', 'pvalue']

plt.figure()
plot = sns.regplot(x="PDB_limit", y="PearsonR", data = data, logx = True)
plot.set_title("Accuracy vs. number of homologs")
plot.set_xlabel("Number of PDB structures used for prediction")
plot.set_ylabel("Pearson correlation of prediction with real B-factors")
plt.savefig('PearsonR_vs_PDB_limit.png')
