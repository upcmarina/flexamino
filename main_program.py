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

# Input fasta file:
querySeq = "rcsb_pdb_4D2I.fasta"  # rcsb pdb header

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
    putative_homologs = blast_my_target(querySeq) # slow execution time

    print("BLAST FINISHED")


# In[ ]:


from input_reader import * 
from pdb_functions import * 
from advanced_models import * 
from profiles import * 
from report import * 
from urllib.request import urlopen

import os
import shutil

querySeq = "rcsb_pdb_4D2I.fasta"  # rcsb pdb header

if __name__ == "__main__":
    
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
    
    profile_predict(querySeq, pdb_paths, "./tmp/alignment.aln", alphafold_model)

    if keep_tmp == False:
        if os.path.exists("./tmp"):
            shutil.rmtree("./tmp")


# In[ ]:


querySeq = "rcsb_pdb_4D2I.fasta"
obtain_uniprot_id(querySeq)


# In[ ]:


import Bio
dir(Bio.PDB.PDBExceptions.PDBConstructionWarning)


# In[ ]:



if __name__ == "__main__":
    
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
    run_clustalw(multifasta)   #output file: ./tmp/alignment.aln

    print("Multiple Sequence Alignment done")

    ######## something is missing

    ### GENERATE ALPHA FOLD MODEL OF QUERY SEQUENCE ###

    alphafold_model = run_alpha_fold(query_uniprot_ID)

    print("Alpha Fold model generated")

    ######## something is missing


# <a style='text-decoration:none;line-height:16px;display:flex;color:#5B5B62;padding:10px;justify-content:end;' href='https://deepnote.com?utm_source=created-in-deepnote-cell&projectId=37b801f1-8b03-4e9d-80a2-bdc088c92f17' target="_blank">
# <img alt='Created in deepnote.com' style='display:inline;max-height:16px;margin:0px;margin-right:7.5px;' src='data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiPz4KPHN2ZyB3aWR0aD0iODBweCIgaGVpZ2h0PSI4MHB4IiB2aWV3Qm94PSIwIDAgODAgODAiIHZlcnNpb249IjEuMSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIiB4bWxuczp4bGluaz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94bGluayI+CiAgICA8IS0tIEdlbmVyYXRvcjogU2tldGNoIDU0LjEgKDc2NDkwKSAtIGh0dHBzOi8vc2tldGNoYXBwLmNvbSAtLT4KICAgIDx0aXRsZT5Hcm91cCAzPC90aXRsZT4KICAgIDxkZXNjPkNyZWF0ZWQgd2l0aCBTa2V0Y2guPC9kZXNjPgogICAgPGcgaWQ9IkxhbmRpbmciIHN0cm9rZT0ibm9uZSIgc3Ryb2tlLXdpZHRoPSIxIiBmaWxsPSJub25lIiBmaWxsLXJ1bGU9ImV2ZW5vZGQiPgogICAgICAgIDxnIGlkPSJBcnRib2FyZCIgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTEyMzUuMDAwMDAwLCAtNzkuMDAwMDAwKSI+CiAgICAgICAgICAgIDxnIGlkPSJHcm91cC0zIiB0cmFuc2Zvcm09InRyYW5zbGF0ZSgxMjM1LjAwMDAwMCwgNzkuMDAwMDAwKSI+CiAgICAgICAgICAgICAgICA8cG9seWdvbiBpZD0iUGF0aC0yMCIgZmlsbD0iIzAyNjVCNCIgcG9pbnRzPSIyLjM3NjIzNzYyIDgwIDM4LjA0NzY2NjcgODAgNTcuODIxNzgyMiA3My44MDU3NTkyIDU3LjgyMTc4MjIgMzIuNzU5MjczOSAzOS4xNDAyMjc4IDMxLjY4MzE2ODMiPjwvcG9seWdvbj4KICAgICAgICAgICAgICAgIDxwYXRoIGQ9Ik0zNS4wMDc3MTgsODAgQzQyLjkwNjIwMDcsNzYuNDU0OTM1OCA0Ny41NjQ5MTY3LDcxLjU0MjI2NzEgNDguOTgzODY2LDY1LjI2MTk5MzkgQzUxLjExMjI4OTksNTUuODQxNTg0MiA0MS42NzcxNzk1LDQ5LjIxMjIyODQgMjUuNjIzOTg0Niw0OS4yMTIyMjg0IEMyNS40ODQ5Mjg5LDQ5LjEyNjg0NDggMjkuODI2MTI5Niw0My4yODM4MjQ4IDM4LjY0NzU4NjksMzEuNjgzMTY4MyBMNzIuODcxMjg3MSwzMi41NTQ0MjUgTDY1LjI4MDk3Myw2Ny42NzYzNDIxIEw1MS4xMTIyODk5LDc3LjM3NjE0NCBMMzUuMDA3NzE4LDgwIFoiIGlkPSJQYXRoLTIyIiBmaWxsPSIjMDAyODY4Ij48L3BhdGg+CiAgICAgICAgICAgICAgICA8cGF0aCBkPSJNMCwzNy43MzA0NDA1IEwyNy4xMTQ1MzcsMC4yNTcxMTE0MzYgQzYyLjM3MTUxMjMsLTEuOTkwNzE3MDEgODAsMTAuNTAwMzkyNyA4MCwzNy43MzA0NDA1IEM4MCw2NC45NjA0ODgyIDY0Ljc3NjUwMzgsNzkuMDUwMzQxNCAzNC4zMjk1MTEzLDgwIEM0Ny4wNTUzNDg5LDc3LjU2NzA4MDggNTMuNDE4MjY3Nyw3MC4zMTM2MTAzIDUzLjQxODI2NzcsNTguMjM5NTg4NSBDNTMuNDE4MjY3Nyw0MC4xMjg1NTU3IDM2LjMwMzk1NDQsMzcuNzMwNDQwNSAyNS4yMjc0MTcsMzcuNzMwNDQwNSBDMTcuODQzMDU4NiwzNy43MzA0NDA1IDkuNDMzOTE5NjYsMzcuNzMwNDQwNSAwLDM3LjczMDQ0MDUgWiIgaWQ9IlBhdGgtMTkiIGZpbGw9IiMzNzkzRUYiPjwvcGF0aD4KICAgICAgICAgICAgPC9nPgogICAgICAgIDwvZz4KICAgIDwvZz4KPC9zdmc+' > </img>
# Created in <span style='font-weight:600;margin-left:4px;'>Deepnote</span></a>
