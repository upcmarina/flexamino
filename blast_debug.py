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
    print("BLAST IS RUNNING...")
    putative_homologs = blast_my_target(querySeq) # slow execution time

    print("BLAST FINISHED")

    pickle.dump(putative_homologs, open("putative_homologs.dat", "wb"))
