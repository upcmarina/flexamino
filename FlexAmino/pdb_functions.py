#!/usr/bin/env python
# coding: utf-8

"""
Package with functions to deal with pdb.
by: Toro, Vallejo, Vega
2022
"""
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio import PDB, SeqIO
from Bio.SeqUtils import IUPACData
import seaborn as sns
import matplotlib.pyplot as plt
import requests as rq
import os

# def searchPDB_crystals(pfam_code, min_res=2.2):
#     """Search PDB crystallographic structures of a given pfam domain below a given resolution threshold."""

#     base_url = "https://search.rcsb.org/rcsbsearch/v1/query?json="
#     pdb_query = '{"query": {"type": "group","logical_operator": "and","nodes": [{"type": "terminal","service": "text","parameters": {"operator": "exact_match","value": "X-RAY DIFFRACTION","attribute": "exptl.method"}},{"type": "terminal","service": "text","parameters":{"operator":"less_or_equal","value": '+ str(min_res) +',"attribute": "rcsb_entry_info.resolution_combined"}},{"type": "terminal","service": "text","parameters": {"operator": "exact_match","value": "'+ pfam_code +'","attribute": "rcsb_polymer_entity_annotation.annotation_id"}}]},"request_options":{"return_all_hits": true},"return_type": "entry"}'
#     get_query = rq.get(base_url+pdb_query)
#     get_query.raise_for_status()
#     response = get_query.json()
#     matches_list = list()  # limitar n results

#     for match in range(0, len(response['result_set'])):
#         matches_list.append(response['result_set'][match]['identifier'])
#     return matches_list


def pdb_download_structure(code): # que descarregui .gz (opcional)
    """Download a specified PDB file using API."""
    url = "https://files.rcsb.org/download/" + code.upper() + ".pdb"
    r = rq.get(url, allow_redirects=True)
    r.raise_for_status()
    PDB_data = r.content.decode("utf-8")

    PDBfile = open("./tmp/"+code.upper()+".pdb", 'wt')
    for line in PDB_data:
        PDBfile.write(line)
    PDBfile.close()


                # afegir generator sara


                # crear nova funcio per cadenes

    ### LINES TO SPLIT THE PDB IN CHAINS
    #parser=PDB.PDBParser()
    #io=PDB.PDBIO()

    #structure = parser.get_structure(out, code.upper()+".pdb")
    #chainnum = 0
    #chainsData = []

    #for chain in structure.get_chains():
    #    chainnum += 1
     #   if chainnum == input_chainnum or chain.get_id() == input_chainnum:
      #      io.set_structure(chain)
       #     io.save("static/data/"+out+".pdb")

    #os.remove("static/data/tempPDB.pdb")

def big_pdb(crystal_dict,limit=50):
    """Download PDB files from a list, with a limit to the number of files."""
    for key in crystal_dict:
        print(key) # check, delete later
        i = 0
        for pdb_id in crystal_dict[key]:
            if i < limit:
                pdb_download_structure(pdb_id) # agafa familiars de la llista i descarregals
                i += 1

def PDB_to_fasta(code_list, inputSeq):
    """Write a FASTA with template sequences from their PDBs and with target FASTA sequence."""
    parser=PDB.PDBParser()
    out_file = open("./tmp/seqs.fasta", 'wt')
    for code in code_list:
        sequence = ""
        structure = parser.get_structure(code[-10:-4], code)
        heteroatom_count= 0
        for residue in structure.get_residues():
            if residue.get_id()[0] == " ": ## exclude heteroatoms
                resname = residue.get_resname()
                one_letter_res = IUPACData.protein_letters_3to1[resname.capitalize()]
                sequence = sequence + one_letter_res
            else:
                heteroatom_count += 1
        print("Excluded heteroatoms:",heteroatom_count ) # change later
        out_file.write(">"+code+"\n")
        out_file.write(sequence+"\n")
    inputSequence = SeqIO.parse(inputSeq, 'fasta')
    for seq in inputSequence:
        out_file.write(str(">"+seq.id+"\n"))
        out_file.write(str(seq.seq+"\n"))
    out_file.close()
    return "./tmp/seqs.fasta"


def pdb_download_chain_xray(code, chain): # que descarregui .gz (opcional)
    """
    Download a specified PDB file using API and extract the specified chain.
    Returns the path to the one-chain pdb file.
    WARNING: IT ONLY DOWNLOADS X-RAY PDBs. If the pdb is not x-ray, returns None
    """
    url = "https://files.rcsb.org/download/" + code.upper() + ".pdb"
    r = rq.get(url, allow_redirects=True)
    r.raise_for_status()
    PDB_data = r.content.decode("utf-8")

    pdb_path = "./tmp/"+code.upper()+".pdb"

    PDBfile = open(pdb_path, 'wt')
    for line in PDB_data:
        PDBfile.write(line)
    PDBfile.close()

    print("Saved pdb " + pdb_path)

    ## parse pdb to extrat chain
    parser=PDBParser()
    try:
        structure = parser.get_structure(code, pdb_path)
    except PDBConstructionWarning as w:
        print("Warning captured:", w)


    pdb_chain_path = "./tmp/"+code.upper()+ "_" + chain + ".pdb"

    ## check if the structure is X-ray and create pdb with only one chain:
    if structure.header["structure_method"] == 'x-ray diffraction':
        chain_structure = structure[0][chain]
        for residue in chain_structure:
            id = residue.id
            if residue.get_id()[0] != " " or  "CA" not in residue:
                chain_structure.detach_child(id)
        io=PDBIO()
        io.set_structure(chain_structure)
        io.save(pdb_chain_path)
        print("Saved pdb "+ pdb_chain_path)
        return pdb_chain_path
    else:
        print(code + " is not x-ray")


def obtain_pdb_list(blast_record):
    """Parses a blast-record object and returns a generator of tuples with the structure (pdb_id, chain) """
    for hit in blast_record.descriptions:
        blast_id = hit.accession
        yield tuple(blast_id.split("_"))
