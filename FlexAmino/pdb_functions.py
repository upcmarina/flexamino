#!/usr/bin/env python
# coding: utf-8

"""
Package with functions to deal with pdbs.
by: Toro, Vallejo, Vega
2022
"""
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio import PDB, SeqIO
from Bio.SeqUtils import IUPACData
import requests as rq
import os
import sys

def PDB_to_fasta(code_list, inputSeq):
    """
    Write a multiFASTA file with template sequences from their PDBs and with target FASTA sequence.

    Arguments:
        code_list:  list of paths to the pdb files whose sequences will be included in the multiFASTA file
        inputSeq:   path to the fasta file containing the query sequence

    Returns:
        A string containing the path to the created multiFASTA file.
    """
    parser=PDB.PDBParser(QUIET=True)
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
        #print("Excluded heteroatoms:",heteroatom_count ) # change later
        out_file.write(">"+code+"\n")
        out_file.write(sequence+"\n")
    inputSequence = SeqIO.parse(inputSeq, 'fasta')
    for seq in inputSequence:
        out_file.write(str(">"+seq.id+"\n"))
        out_file.write(str(seq.seq+"\n"))
    out_file.close()
    return "./tmp/seqs.fasta"


def pdb_download_chain_xray(code, chain, verbose): # que descarregui .gz (opcional)
    """
    Download a specified PDB file using API and extract the specified chain.

    Arguments:
        code:   a pdb code
        chain:  a valid chain ID
        verbose: if True, it will print the messages on the screen, if False it will not.

    Returns:
        The path to the one-chain pdb file.
        WARNING: IT ONLY DOWNLOADS X-RAY PDBs. If the pdb is not x-ray, returns None.
    """
    pdb_path = "./tmp/"+code.upper()+".pdb"

    if os.path.exists(pdb_path):
        if verbose: print(pdb_path + " already exists.", file=sys.stderr, flush=True)
    else:
        url = "https://files.rcsb.org/download/" + code.upper() + ".pdb"

        try:
            r = rq.get(url, allow_redirects=True)
            r.raise_for_status()
        except rq.exceptions.HTTPError:
            if verbose: print("PDB " + code.upper() + " not found", file=sys.stderr, flush=True)
            return None

        PDB_data = r.content.decode("utf-8")

        PDBfile = open(pdb_path, 'wt')
        for line in PDB_data:
            PDBfile.write(line)
        PDBfile.close()

        if verbose: print("Saved pdb " + pdb_path, file=sys.stderr, flush=True)


    ## parse pdb to extract chain
    parser=PDBParser(QUIET=True)
    structure = parser.get_structure(code, pdb_path)


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
        if verbose: print("Saved pdb "+ pdb_chain_path, file=sys.stderr, flush=True)
        return pdb_chain_path
    else:
        if verbose: print(code + " is not x-ray", file=sys.stderr, flush=True)


def obtain_pdb_list(blast_record):
    """
    Parses the hits of a blast-record object and returns a generator of tuples with the structure (pdb_id, chain) 

    Arguments:
        blast_record:   blast-record object

    Returns:
        Generator of tuples with the structure (pdb_id, chain_id)

    """
    for hit in blast_record.descriptions:
        blast_id = hit.accession
        yield tuple(blast_id.split("_"))
