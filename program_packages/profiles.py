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
import prody as pd
import requests as rq
import os
from pathlib import Path


def profile_generator(pdb_list):
    """
    Given a list of .pdb files as an input, 
    for each file compute the average of the b factor values for each amino acid
    """

    ###EXTRACT CA BFACTORS FROM EACH FILE BY AMINOACID
    bfactors = {}

    for pdb in pdb_list:
        
        file_path = "./tmp" +"/"+ pdb
        open_pdb = open(file_path)

        p = PDB.PDBParser()
        pdbname = pdb[:-4] ## to remove ".pdb"
        structure = p.get_structure(pdbname, open_pdb)
        bfactors[pdbname] = {}
  
        for residue in structure.get_residues():
  
            if residue.get_id()[0] == " ": ## exclude heteroatoms
                resname = residue.get_resname()
                shortresname = IUPACData.protein_letters_3to1[resname.capitalize()] 
                CA_bfactor = residue["CA"].get_bfactor()
                bfactors[pdbname].setdefault(shortresname, []).append(CA_bfactor)

    ##### NORMALIZE BFACTORS BEFORE COMPUTING THE MEANS???

    ## COMPUTE AVERAGE FOR EACH AA FOR EACH FILE
    average_bfactors = {}

    for pdb in bfactors.keys():
  
        for aa, bfactorlist in bfactors[pdb].items():
            mean = sum(bfactorlist)/len(bfactorlist)
            average_bfactors.setdefault(aa, []).append(mean)

    ## NOW compute the mean between files ????
    bfactors_profile = {}

    for aa, mean_bfactors in average_bfactors.items():
        bfactors_profile[aa] = sum(mean_bfactors)/len(mean_bfactors)

    return bfactors_profile   # what, a python obj like dict , later take that dict and pass it to file


def profile_predict(inputSeq, family_profile):
    """Predict b-factor profile for a given aminoacid sequence.
    
    Takes an input an aminoacid sequence and a previously generated beta-factor profile
    to predict beta factors for that sequence.

    Keyword arguments:
    inputSeq -- the input aminoacid sequence (as a string) for which the b-factor prediction will be made.
    family_profile -- the family beta-factor profile (as a dictionary) which is applied to the input sequence.
    """
  
    bf_profile = list()
   
    for seq_record in SeqIO.parse(inputSeq, "fasta"):
        input_seq =  str(seq_record.seq)

    for AA in input_seq:
        bf_profile.append(family_profile[AA])

    return (input_seq, bf_profile)
