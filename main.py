""" 
Main execution program. 
by: Toro, Vallejo, Vega
2022
"""

from input_reader import * 
from pdb_functions import * 
from profiles import * 
from report import * 

import os
import shutil

querySeq = "example.fa"

if __name__ == "__main__":

    ### INPUT READER ###

    significant_families = getPfamilies(querySeq)
    print(significant_families) # check, delete later


    ### PDB FUNCTIONS ###

    crystal_dict = {} # create empty dictionary
    for key in significant_families:
        crystal_dict[key] = searchPDB_crystals(key)
    print(crystal_dict) # check, delete later

    if os.path.exists("./tmp"):
        shutil.rmtree("./tmp")

    os.mkdir("./tmp") # create temporary directory to store .pdb files

    big_pdb(crystal_dict,limit=3)


    ### PROFILES ###

    # Train:
    pdb_list = os.listdir('./tmp')
    print(pdb_list) # check, delete later
    trained_profiles = profile_generator(pdb_list)
    print(trained_profiles) # check, delete later
    
    # Scoring:
    my_profile = profile_predict(querySeq,trained_profiles)
    print(my_profile) # check, delete later

    ### OUTPUT ####
    
    # Obtain parseable .txt file with flexibility scores:
    prediction_write(my_profile,"my_profile_file.txt") 
    
    # Obtain plots:
    plot_profile(my_profile,querySeq)
