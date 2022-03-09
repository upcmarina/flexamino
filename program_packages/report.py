#!/usr/bin/env python
# coding: utf-8

""" 
Package with functions to generate a parseable file with the results and plots. 
by: Toro, Vallejo, Vega
2022
"""

from Bio import PDB, SeqIO
from Bio.SeqUtils import IUPACData
import seaborn as sns
import matplotlib.pyplot as plt
import prody as pd
import requests as rq
import os
import plotly.express as px



def prediction_write(prediction, output_filename, pdb=None, chainID=None):
    """
    Write predicted beta-factors to a file.
    Takes a predicted beta-factor profile and writes it to a file.
    If a PDB file is passed, it writes the predicted beta-factors to the alpha carbon's beta-factor column;
    if no PDB is given, writes the output in a tabular format.
    Keyword arguments:
    prediction -- the tuple containing the aminoacid sequence and the predicted beta-factors
    output_filename -- the name with which to write the output file
    pdb -- the name of the pdb to modify
    """
 
    if pdb is None:
        with open(output_filename, "wt") as fd:
            fd.write(">"+output_filename+"\n")
            for position in range(0, len(prediction[0])):
                fd.write(prediction[0][position]+"\t"+str(prediction[1][position])+"\n")
    #else: # CAL PODER FICAR ELS BFACTORS A UN PDB?
     #   structure = p.get_structure("structure", pdb)
      #  for chain in structure.get_chains():
       #     if chain == chainID:
        #        position = 0
         #       for residue in structure.get_residues():
          #          #resname = residue.get_resname()
           #         for atom in residue.get_atoms():
            #            if atom.get_id() == "CA":
             #               atom.set_bfactor(prediction[1][position])
              #              position += 1

    #return structure


def plot_profile(prediction, seqID):
    """Plot the predicted b-factor scores for each residue position."""

    AA_pos = [num for num in range(0, len(prediction[0]))]
    plt.figure()
    plot = sns.lineplot(x = AA_pos, y = prediction[1])
    plot.set_title("Predicted B-factor profile for sequence " + seqID)
    plot.set_xlabel("Residue number")
    plot.set_ylabel("Normalized B-factor")

def interactive_plot_profile(prediction, seqID):
    AA_pos = [num for num in range(0, len(prediction[0]))]
    fig = px.line(data=AA_pos,x="AA_pos", y="prediction[1]")
    fig.show()
