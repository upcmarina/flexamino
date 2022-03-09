#!/usr/bin/env python
# coding: utf-8

""" 
Package with functions to deal with pdb. 
by: Toro, Vallejo, Vega
2022
"""

from Bio import PDB, SeqIO
from Bio.SeqUtils import IUPACData
import seaborn as sns
import matplotlib.pyplot as plt
import prody as pd
import requests as rq
import os

def searchPDB_crystals(pfam_code, min_res=2.2):
    """Search PDB crystallographic structures of a given pfam domain below a given resolution threshold."""
 
    base_url = "https://search.rcsb.org/rcsbsearch/v1/query?json="
    pdb_query = '{"query": {"type": "group","logical_operator": "and","nodes": [{"type": "terminal","service": "text","parameters": {"operator": "exact_match","value": "X-RAY DIFFRACTION","attribute": "exptl.method"}},{"type": "terminal","service": "text","parameters":{"operator":"less_or_equal","value": '+ str(min_res) +',"attribute": "rcsb_entry_info.resolution_combined"}},{"type": "terminal","service": "text","parameters": {"operator": "exact_match","value": "'+ pfam_code +'","attribute": "rcsb_polymer_entity_annotation.annotation_id"}}]},"request_options":{"return_all_hits": true},"return_type": "entry"}'
    get_query = rq.get(base_url+pdb_query)
    get_query.raise_for_status()
    response = get_query.json()
    matches_list = list()  # limitar n results

    for match in range(0, len(response['result_set'])):
        matches_list.append(response['result_set'][match]['identifier'])
    return matches_list
