#!/usr/bin/env python3
# coding: utf-8

"""
setup.py
Authors : Toro; Vallejo; Vega
2022
"""

from setuptools import setup
from os import path

setup(name='FlexAmino',
  version='1.0',
  description='FlexAmino, Tool to obtain Flexibility Scores for each aminoacid in a given sequence',
  long_description=open('README.md').read(),
  author='Eric Toro Delgado, Marina Vallejo Vallés, Sara Vega Abellaneda',
  author_email='eric.toro01@estudiant.upf.edu, marina.vallejo01@estudiant.upf.edu, sara.vega02@estudiant.upf.edu',
  url='https://github.com/upcmarina/structural_biology',
  install_requires = ['Biopython'],
  packages=['FlexAmino'],
  scripts=['FlexAmino/flexamino.py', 'FlexAmino/args.py', 'FlexAmino/pdb_functions.py', 'FlexAmino/profiles_bfactor.py','FlexAmino/work_with_input_sequence.py', 'FlexAmino/results_to_file.py'],
  keywords = 'bioinformatics protein_flexibility b_factors'
)
