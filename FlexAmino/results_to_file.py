#!/usr/bin/env python
# coding: utf-8

"""
Package with functions to generate a parseable file with the results and plots.
by: Toro, Vallejo, Vega
2022
"""

import seaborn as sns
import matplotlib.pyplot as plt


def prediction_write(prediction, output_filename, seqID, winsize):
    """
    Takes a predicted beta-factor profile and writes it to a file.

    Arguments:
        prediction: the tuple containing the aminoacid sequence and the predicted beta-factors
        output_filename:    prefix of the output filename
        seqID:  ID of the query sequence
        winsize:    size of the sliding window applied in the smoothing process

    Returns:
        None
    """

    with open(output_filename + ".txt", "wt") as fd:
        fd.write(">"+seqID+"\n")
        if winsize == 1:
            for position in range(0, len(prediction[0])):
                fd.write(str(position+1)+"\t"+prediction[0][position]+
                "\t"+str(prediction[1][position])+
                "\t"+str(prediction[2][position])+"\n")
        else:
            for position in range(0, len(prediction[0])):
                fd.write(str(position+1)+"\t"+prediction[0][position]+
                "\t"+str(prediction[1][position])+
                "\t"+str(prediction[2][position])+
                "\t"+str(prediction[3][position])+"\n")


def plot_profile(prediction, seqID, output_filename, winsize):
    """
    Plot the predicted b-factor scores for each residue position.

    Arguments:
        prediction: the tuple containing the aminoacid sequence and the predicted beta-factors
        output_filename:    prefix of the output filename
        seqID:  ID of the query sequence
        winsize:    size of the sliding window applied in the smoothing process

    Returns:
        None
    """

    AA_pos = [num for num in range(0, len(prediction[0]))]

    plt.figure()
    plot = sns.lineplot(x = AA_pos, y = prediction[3])
    if winsize > 1:
        plot.set_title("Predicted B-factor profile for sequence " + seqID + "\n(window size " + str(winsize) + ")")
    else:
        plot.set_title("Predicted B-factor profile for sequence " + seqID)
    plot.set_xlabel("Residue number")
    plot.set_ylabel("Normalized B-factor")

    if output_filename is not None:
        plt.savefig(output_filename + '.png')
