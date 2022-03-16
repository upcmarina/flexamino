#!/usr/bin/env python
# coding: utf-8

# # Main program, execution

# In[ ]:


"""
Main execution program.
by: Toro, Vallejo, Vega
2022
"""

from input_reader import *
from pdb_functions import *
from advanced_models import *
from profiles import *
from report import *
from urllib.request import urlopen

import os
import shutil
import pickle

# Input fasta file:
querySeq = "Q16613.fasta"  # rcsb pdb header
output_filename = "output.out"
# Do you want to keep the tmp directory?
keep_tmp = True

if __name__ == "__main__":

    # Create temporary directory to store .pdb files
    # if directory exists, delete it:
    if os.path.exists("./tmp"):
        shutil.rmtree("./tmp")
    os.mkdir("./tmp")

    ### INPUT READER ###

    # obtain UNIPROT id:
    query_uniprot_ID = obtain_uniprot_id(querySeq)

    # run BLAST:
    # print("BLAST IS RUNNING...")
    # putative_homologs = blast_my_target(querySeq) # slow execution time
    #
    # print("BLAST FINISHED")

    putative_homologs = pickle.load(open("putative_homologs.dat", "rb"))


    ### Parse blast output AND download pdbs in tmp directory AND save filepaths in a list
    pdb_paths = [] #list with the paths to the downloaded pdb files
    pdb_limit = 3 #user can change this n,
    for pdb_hit in obtain_pdb_list(putative_homologs):
        if len(pdb_paths) >= pdb_limit:
            break
        pdb_id = pdb_hit[0]
        chain = pdb_hit[1]
        path = pdb_download_chain_xray(pdb_id, chain)
        if path != None:  ## if the pdb is not x-ray, the function pdb_download_chain returns None
            pdb_paths.append(path)

    # obtain multifasta file:
    multifasta = PDB_to_fasta(pdb_paths, querySeq)

    print("Multifasta file obtained")

    # Multiple Sequence Alignment (MSA)
    msa = run_clustalw(multifasta)

    print("Multiple Sequence Alignment done")

    ######## something is missing

    ### GENERATE ALPHA FOLD MODEL OF QUERY SEQUENCE ###

    alphafold_model = run_alpha_fold(query_uniprot_ID)

    print("Alpha Fold model generated")

    ######## something is missing

    prediction = profile_predict(querySeq, pdb_paths, msa, alphafold_model)

    prediction_write(prediction, output_filename)

    if keep_tmp == False:
        if os.path.exists("./tmp"):
            shutil.rmtree("./tmp")
