#!/usr/bin/env python
# coding: utf-8

""" 
Package with functions to deal with input fasta files. 
by: Toro, Vallejo, Vega
2022
"""

from Bio import PDB, SeqIO
from Bio.SeqUtils import IUPACData
import prody as pd
import requests as rq

def getPfamilies(querySeq):
    """Search the input sequence in Pfam to get the matching protein families."""
    for seq_record in SeqIO.parse(querySeq, "fasta"):
        query_seq =  str(seq_record.seq)

    pfamilies = pd.searchPfam(query_seq) # llegeix una string
    families = {}
    for domain in pfamilies.keys():
        families[domain] = {'align_start': pfamilies[domain]['locations']['ali_start'], 'align_end': pfamilies[domain]['locations']['ali_end']}

    return families
