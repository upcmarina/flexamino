#!/usr/bin/env python3
# coding: utf-8

# # Main program, execution

# In[ ]:


"""
Main execution program.
by: Toro, Vallejo, Vega
2022
"""
import args

from input_reader import *
from pdb_functions import *
from advanced_models import *
from profiles import *
from report import *

from urllib.request import urlopen

import os
import shutil
import pickle
import sys
import time

# Proteins to test: P06401, Q9P7Q4, Q9VVG4, P38401, P11433, Q9Y223, P65206

if __name__ == "__main__":
    
    start = time.time()

    # Get command-line arguments
    options = args.getArgs()
    querySeq = options.infile
    if querySeq == None:
        raise SystemExit('Please provide an input file')
    output_filename = options.outfile
    keep_tmp = options.keep_tmp
    verbose = options.verbose
    rescue = options.rescue
    pdb_limit = options.pdb_limit
    
    
    print("### FlexAmino HAS STARTED", file=sys.stderr, flush=True)

    # Create temporary directory to store .pdb files
    # if directory exists, delete it:
    try:
        os.mkdir("./tmp")
    except FileExistsError:
        pass

    ### INPUT READER ###

    # obtain UNIPROT id:
    query_uniprot_ID = obtain_uniprot_id(querySeq)

    if rescue == False:
    # run BLAST:
        if verbose: print("### BLAST is running...", file=sys.stderr, flush=True) 
        putative_homologs = blast_my_target(querySeq) # slow execution time
        with open("./tmp/putative_homologs.dat", 'wb') as blast_file:
            pickle.dump(putative_homologs, blast_file)
    #
        if verbose: print("### BLAST finished", file=sys.stderr, flush=True)

    else:
        if verbose: print("### Recovering BLAST results...", file=sys.stderr, flush=True)
        putative_homologs = pickle.load(open("./tmp/putative_homologs.dat", "rb"))


    ### Parse blast output AND download pdbs in tmp directory AND save filepaths in a list
    pdb_paths = [] #list with the paths to the downloaded pdb files
    
    for pdb_hit in obtain_pdb_list(putative_homologs):
        if len(pdb_paths) >= pdb_limit:
            break
        pdb_id = pdb_hit[0]
        chain = pdb_hit[1]
        path = pdb_download_chain_xray(pdb_id, chain, verbose)
        if path != None:  ## if the pdb is not x-ray, the function pdb_download_chain returns None
            pdb_paths.append(path)

    # obtain multifasta file:
    multifasta = PDB_to_fasta(pdb_paths, querySeq)

    if verbose: print("### Multifasta file obtained", file=sys.stderr, flush=True)

    # Multiple Sequence Alignment (MSA)
    msa = run_clustalw(multifasta)

    if verbose: print("### Multiple Sequence Alignment done", file=sys.stderr, flush=True)

    ######## something is missing

    ### GENERATE ALPHA FOLD MODEL OF QUERY SEQUENCE ###

    alphafold_model = run_alpha_fold(query_uniprot_ID)

    if verbose: print("### Alpha Fold model generated", file=sys.stderr, flush=True)

    ######## something is missing

    prediction = profile_predict(querySeq, pdb_paths, msa, alphafold_model)

    prediction_write(prediction, output_filename, query_uniprot_ID)
    
    if verbose: print("### Output generated", file=sys.stderr, flush=True)

    plot_profile(prediction, query_uniprot_ID, output_filename)

    if keep_tmp == False:
        if os.path.exists("./tmp"):
            shutil.rmtree("./tmp")

        if verbose: print("### tmp directory deleted", file=sys.stderr, flush=True)

    end = time.time()

    print("### FlexAmino has finished. Total execution time: ", end-start, " s", file=sys.stderr, flush=True)
