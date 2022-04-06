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

from work_with_input_sequence import *
from pdb_functions import *
from profiles_bfactor import *
from results_to_file import *

from urllib.request import urlopen

import os
import shutil
import pickle
import sys
import time

if __name__ == "__main__":
    
    start = time.time()

    print("### FlexAmino HAS STARTED", file=sys.stderr, flush=True)

    # Get command-line arguments:
    options = args.getArgs()
    querySeq = options.infile
    if querySeq == None:
        raise SystemExit('Please provide an input file.')
    output_filename = options.outfile
    if output_filename == None:
        raise SystemExit('Please provide an output prefix.')
    keep_tmp = options.keep_tmp
    verbose = options.verbose
    rescue = options.rescue
    pdb_limit = options.pdb_limit
    winsize = options.winsize
 
    # Create temporary directory to store .pdb files
    # if directory exists, delete it:
    try:
        os.mkdir("./tmp")
    except FileExistsError:
        pass

    # Obtain UNIPROT id from the Query Sequence:
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

    # Obtain multifasta file:
    multifasta = PDB_to_fasta(pdb_paths, querySeq)

    if verbose: print("### Multifasta file obtained", file=sys.stderr, flush=True)

    # Multiple Sequence Alignment (MSA):
    msa = run_clustalw(multifasta)

    if verbose: print("### Multiple Sequence Alignment done", file=sys.stderr, flush=True)

    # Generate Alpha Fold model of the query sequence
    alphafold_model = run_alpha_fold(query_uniprot_ID)

    if verbose: print("### Alpha Fold model downloaded", file=sys.stderr, flush=True)

    prediction = profile_predict(querySeq, pdb_paths, msa, alphafold_model, winsize)

    prediction_write(prediction, output_filename, query_uniprot_ID, winsize)
    
    if verbose: print("### Output generated", file=sys.stderr, flush=True)

    plot_profile(prediction, query_uniprot_ID, output_filename)

    if keep_tmp == False:
        if os.path.exists("./tmp"):
            shutil.rmtree("./tmp")

        if verbose: print("### tmp directory deleted", file=sys.stderr, flush=True)

    end = time.time()

    print("### FlexAmino has finished. Total execution time: ", end-start, " s", file=sys.stderr, flush=True)
