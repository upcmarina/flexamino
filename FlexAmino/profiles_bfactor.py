#!/usr/bin/env python
# coding: utf-8

"""
Package with functions to generate and apply b-factor profiles.
by: Toro, Vallejo, Vega
2022
"""

from Bio import PDB, SeqIO, AlignIO
import numpy as np
import warnings


def normalize_bfactors(pdb_list):
    """
    Normalize c-alfa b-factors of a PDB structure by computing their Z-score.

    Arguments:
        pdb_list:   list of paths of the PDB files to be normalized
    Returns:
        List of objects of the class Structure that corresponds to the parsed input PDBs with the normalized b-factors.
    """
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


def profile_predict(target_name, templates_list, fasta_align, alphaFold_path, winsize):
    """
    Predict the b-factor profile of a target protein based on a structural alignment.

    Arguments:
        target_name:    path to the query sequence fasta file
        templates_list: list of paths of the PDB files to be used as templates
        fasta_align:    path to the MSA file in clustal format of the templates and the query sequences. The query sequence must be in the last position
        alphaFold_path: path to the PDB structure of the query sequence
        winsize:        sliding window size to be applied in the smoothing of the predicted b-factors

    Returns:
        A tuple containing the aminoacid sequence, the predicted beta-factors, the standard deviation of the beta-factors and the smoothed b-factors.
    """
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


        position = 0
        for residue in structural_aln.get_maps()[0].values():
            if residue is None:
                CA_bfactors[element, position] = np.nan
            else:
                CA_bfactors[element, position] = residue['CA'].get_bfactor()
            position += 1


        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mean_bfactors = np.nanmean(CA_bfactors, axis = 0)
            std_bfactors = np.nanstd(CA_bfactors, axis = 0)

            smooth_mean = np.empty(len(mean_bfactors))
            smooth_std = np.empty(len(std_bfactors))

            for position in range(0, len(mean_bfactors)):
                if winsize % 2 != 0:
                    lowerbound = int(position - (winsize - 1)/2)
                else:
                    lowerbound = int(position - winsize/2 + 1)

                upperbound = int(lowerbound + winsize)

                if (lowerbound >= 0 and upperbound <= len(mean_bfactors) and not np.isnan(mean_bfactors[position])):
                    smooth_mean[position] = np.nanmean(mean_bfactors[lowerbound:upperbound])
                else:
                    smooth_mean[position] = np.nan


    return (str(target_seq.seq), mean_bfactors.tolist(), std_bfactors.tolist(), smooth_mean.tolist())
