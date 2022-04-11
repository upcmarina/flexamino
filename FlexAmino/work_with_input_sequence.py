#!/usr/bin/env python
# coding: utf-8

"""
Package with functions to deal with the input fasta file, perform BLAST, CLUSTALW and generate an Alpha Fold model from it.
by: Toro, Vallejo, Vega
2022
"""

from Bio import PDB, SeqIO
from Bio.SeqUtils import IUPACData
from Bio.Blast.NCBIWWW import qblast
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline

import os
import requests as rq
import urllib.parse
import urllib.request
from urllib.request import urlopen
import shutil

def blast_my_target(inputSeq):
    """
    BLAST
    """
    fasta_file = open(inputSeq, "r")
    seq = fasta_file.read()
    fasta_file.close()

    blast_result = qblast(url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi', \
        format_object='Alignment', \
        program="blastp", \
        database="pdb", \
        sequence=seq, \
        format_type="XML")

    blast_record = NCBIXML.read(blast_result)

    return blast_record

def run_clustalw(multifasta):
    """
    Performs a Multiple Sequence Alignment using ClustalW.
    Accepts a multifasta file as input.
    Returns a .aln file
    """
    file_path = os.path.dirname(os.path.abspath(__file__))
    os.system(file_path + './bin/clustalw2 ' + multifasta + ' -OUTORDER=INPUT > ./tmp/clustalw.log')
    return "./tmp/seqs.aln"


def obtain_uniprot_id(inputSeq):
    """
    Obtain uniprot id of our query sequence
    """
    with open(inputSeq) as my_file:

        records = SeqIO.parse(my_file, 'fasta')

        for record in records:

            # Get id from a uniprot file:
            if record.id[:2]=="sp":
                return record.id[3:9]

            # Get id from pdb file and convert to uniprot:
            else:
                url = 'https://www.uniprot.org/uploadlists/'
                params = {
                    'from': 'PDB_ID',
                    'to': 'ACC',
                    'format': 'tab',
                    'query': record.id[0:4]
                    }

                data = urllib.parse.urlencode(params)
                data = data.encode('utf-8')
                req = urllib.request.Request(url, data)

                with urllib.request.urlopen(req) as f:
                    response = f.read()
                    code=response.decode('utf-8')
                    return(code[-7:-1])


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
