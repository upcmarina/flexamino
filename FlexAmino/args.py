import argparse
import sys
import os
import Bio.Seq


def getArgs():
    parser = argparse.ArgumentParser(description = "Open protein FASTA files")
    parser.add_argument('-i', '--input',
                            dest = 'infile', action = 'store',
                            default = None,
                            help = 'An input FASTA file containing protein sequences',
                            type = str)
    parser.add_argument('-o', '--output',
                            dest = 'outfile', action = 'store',
                            default = None,
                            help = 'Set a prefix to generate the different output files.',
                            type = str)
    parser.add_argument('-v', '--verbose',
                            dest = 'verbose', action = 'store_true',
                            default = False,
                            help = 'Print the progression of the program execution to the terminal (Standard Error).')
    parser.add_argument('-t', '--tmp',
                            dest = 'keep_tmp', action = 'store_true', 
                            default = False,
                            help = 'Keep the tmp directory when the program finishes.')
    parser.add_argument('-r', '--rescue',
                            dest = 'rescue', action = 'store_true',
                            default = False,
                            help = 'Recover a computation from a BLAST result to avoid running BLAST again if the program fails.')  
    parser.add_argument('-p', '--pdb_cutoff',
                            dest = 'pdb_limit', action = 'store',
                            default = 10,
                            help = 'Set a maximum number of pdb structures to use for the computation.')  
    parser.add_argument('-w', '--winsize',
                            dest = 'winsize', action = 'store',
                            default = 1,
                            help = 'Set a sliding window for smoothing the results.')  

    return parser.parse_args()









if __name__ == "__main__":
    pass