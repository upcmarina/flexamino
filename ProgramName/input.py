import argparse
import sys
import os
import Bio.Seq


def getArgs():
    parser = argparse.ArgumentParser(description = "Open protein FASTA files")
    parser.add_argument('-i', '--input',
                            dest = 'infile', action = 'store',
                            default = './'
                            help = 'An input FASTA file containing protein sequences (or a directory?????)',
                            type = str)
    parser.add_argument('-o', '--output',
                            dest = 'oufile', action = 'store',
                            help = 'Name to use to generate the different output files',
                            type = str)
    parser.add_argument('-v', '--verbose',
                            dest = 'verbose', action = 'store_true',
                            default = False,
                            help = 'Print the progression of the program to the terminal (Standard Error)')
    return parser.parse_args()









if __name__ == "__main__":
