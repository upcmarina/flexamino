#!/usr/bin/env python
# coding: utf-8

""" 
Package with functions to deal with Alpha fold. 
by: Toro, Vallejo, Vega
2022
"""

from input_reader import * 
from pdb_functions import * 
from profiles import * 
from report import * 
from urllib.request import urlopen

import os
import shutil

def run_alpha_fold(query_uniprot):
    """
    Generate an alpha fold model of the query sequence.
    """
    url = "https://alphafold.ebi.ac.uk/files/AF-"+query_uniprot+"-F1-model_v2.pdb"
    filehandle = urlopen(url)
    pdbfile = open("./tmp/"+ query_uniprot +"_alphaFold.pdb", 'wb')
    
    for line in filehandle:
        pdbfile.write(line)
    pdbfile.close()

    return "./tmp/"+ query_uniprot +"_alphaFold.pdb"