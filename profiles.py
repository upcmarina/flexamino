#!/usr/bin/env python
# coding: utf-8

""" 
Package with functions to generate and apply b-factor profiles. 
by: Toro, Vallejo, Vega
2022
"""

from Bio import PDB, SeqIO
from Bio.SeqUtils import IUPACData
import seaborn as sns
import matplotlib.pyplot as plt
import requests as rq
import os
from pathlib import Path


def normalize_bfactors(pdb_list):
    """Normalize c-alfa b-factors of a PDB structure by computing their Z-score."""
    normalized_structures = list()
    for pdb_name in pdb_list:
        structure = parser.get_structure(pdb_name)
        bfactors = list()
        for residue in structure.get_residues():
            bfactors.append(residue['CA'].get_bfactor())
        bf_mean = np.mean(bfactors)
        bf_sd = np.std(bfactors)
        for residue in structure.get_residues():
            residue['CA'].set_bfactor((residue['CA'].get_bfactor() - bf_mean)/bf_sd)
        normalized_structures.append(structure)
    return normalized_structures


def profile_predict(target_name, templates_list, fasta_align, alphaFold_path):
    """Predict the b-factor profile of an target protein based on a structural alignment."""
    parser=PDB.PDBParser()

    target_seq = SeqIO.read(target_name, "fasta")

    CA_bfactors = np.empty((0, len(target_seq)))

    target_structure = parser.get_structure(target_name[-10:-4], alphaFold_path)

    for element in range(0, len(templates_list)):
        template_structure = parser.get_structure(templates_list[element][-10:-4], templates_list[element])
        template_name = template_structure.get_id()
        structural_aln = PDB.StructureAlignment(fasta_align, target_structure, template_structure, len(templates_list), element)
        
        bfactor_list = np.empty(0)
        
        position = 0
        for residue in structural_aln.get_maps()[0].values():
            if residue is None:
                bfactor_list = np.append(bfactor_list, [np.nan], axis = 0)
                #print(residue)
            else:
                bfactor_list = np.append(bfactor_list, [residue['CA'].get_bfactor()], axis = 0)
                #print(residue)
            position += 1

        CA_bfactors = np.append(CA_bfactors, [bfactor_list], axis = 0)


    mean_bfactors = np.nanmean(CA_bfactors, axis=0)

    return (str(target_seq.seq), mean_bfactors.tolist())