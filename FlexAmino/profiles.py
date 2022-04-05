#!/usr/bin/env python
# coding: utf-8

"""
Package with functions to generate and apply b-factor profiles.
by: Toro, Vallejo, Vega
2022
"""

from Bio import PDB, SeqIO, AlignIO
from Bio.SeqUtils import IUPACData
import seaborn as sns
import matplotlib.pyplot as plt
import requests as rq
import os
from pathlib import Path
import numpy as np
import warnings


def normalize_bfactors(pdb_list):
    """Normalize c-alfa b-factors of a PDB structure by computing their Z-score."""
    normalized_structures = list()
    parser=PDB.PDBParser()
    for pdb_name in pdb_list:
        structure = parser.get_structure(pdb_name[-10:-4], pdb_name)
        bfactors = list()
        for residue in structure.get_residues():
            if "CA" in residue:
                bfactors.append(residue['CA'].get_bfactor())
        bf_mean = np.mean(bfactors)
        bf_sd = np.std(bfactors)
        for residue in structure.get_residues():
            if "CA" in residue:
                residue['CA'].set_bfactor((residue['CA'].get_bfactor() - bf_mean)/bf_sd)
        normalized_structures.append(structure)
    return normalized_structures


def profile_predict(target_name, templates_list, fasta_align, alphaFold_path):
    """Predict the b-factor profile of an target protein based on a structural alignment."""
    parser=PDB.PDBParser()

    target_seq = SeqIO.read(target_name, "fasta")

    CA_bfactors = np.empty((len(templates_list), len(target_seq)))

    target_structure = parser.get_structure(target_name[-10:-4], alphaFold_path)
    alignment = AlignIO.read(fasta_align, "clustal")
    norm_templates = normalize_bfactors(templates_list) ## array of structure objects
    for element in range(0, len(templates_list)):
        template_structure = norm_templates[element]
        template_name = template_structure.get_id()
        structural_aln = PDB.StructureAlignment(alignment, target_structure, template_structure, len(templates_list), element)

        #bfactor_list = np.empty(0)

        position = 0
        for residue in structural_aln.get_maps()[0].values():
            if residue is None:
                CA_bfactors[element, position] = np.nan
                #bfactor_list = np.append(bfactor_list, [np.nan], axis = 0)
                #print(residue)
            else:
                CA_bfactors[element, position] = residue['CA'].get_bfactor()
                #bfactor_list = np.append(bfactor_list, [residue['CA'].get_bfactor()], axis = 0)
                #print(residue)
            position += 1

        #CA_bfactors = np.append(CA_bfactors, [bfactor_list], axis = 0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mean_bfactors = np.nanmean(CA_bfactors, axis = 0)
            std_bfactors = np.nanstd(CA_bfactors, axis = 0)

    return (str(target_seq.seq), mean_bfactors.tolist(), std_bfactors.tolist())
