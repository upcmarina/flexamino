**PROJECT INFORMATION**
===================================

# GOAL

Development of flexibility score for proteins based on STRUCTURE. MSA is extra information to take into account and see if flexibility score has sense (in a position without variation in MSA, we can't have a large flexibility score).

Structures from PDB or AlphaFold.

# INPUT

Protein sequence or protein family.

Important, our approach must be general, valid for any type of family and not specific.

# OUTPUT

Flexibility score for each amino acid in sequence. But side chain? or backbone? or both? --> preguntar Baldo
'
We must obtain:

* Parseable text output

* Graphical

# POSSIBLE STEPS

* Multiple Sequence Alignment **MSA**

* Structural Alignment

* Develop flexibility score analyzing alignments

* Graphical representation of flexibility scores (idea: fer-ho com PROSA)

* Test 4 families of proteins (2 proposed by the professors, 2 selected by us)

# OTHER

Indep evaluation for PYT and SBI. Theoretical part --> Baldo

Develop pipeline/ workflow, don't import  a library that calculates the flex scores. Use libraries like biopython, but not packages for flexibility.

To read MSA we can use libraries for that.

Program structure! very important, file for classes, DIVIDE CODE IN MODULES!! 


by analyzing this alignments develop a flex score and then represent results --> we'll see libraries for this scores

Date for project --> near exam day

# IDEAS FLEXIBILITY SCORE

* Use RMSD

* B-factor-derived flexibility score (Smith et al. 2003, Improved amino acid flexibility parameters): from experimental B-facotors of crystallography  
distributions of b-factors for each AA are derived, and the flexiblity score assigned to each AA residue is the location parameter fo the extreme  
valued dsitribution that the b-factors follow for each residue type. We could use this to compute fliexiblity in a sliding-window manner.

* Molecular Dynamics

* Structural alphabets (Dong et al. 2016; Characterization and Prediction of Protein FlexibilityBased on Structural Alphabets):  
Pros:
    - Seems more or less simple to do
    - Does not use B-factor, which they say has experimental noise (depends on resolution of the structure, refinement conditions of the crystal, etc.)  
Cons:
    - From the ROC curves of the article I wouldn't say the performance is great.

# PROGRAM DRAFT
=========================
## Prediction part:

1. Read FASTA sequence(s): this can be done with the following example lines from Python class:
```
from Bio import SeqIO
for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
print(seq_record.id)
print(repr(seq_record.seq))
print(len(seq_record))
```
This makes a generator object for the sequences, like the fasta iterator we did in the exercises.

2. Get the family of the protein: this may be a little bit tricky because one protein can match to multiple families I think? Also I don't know if there is some way to do it with Python modules or we could directly call Pfam/HMMER web services via API (by could I mean if we are allowed to do it, I imagine they will have APIs). Some ideas are:
- Use Pfam website to search the sequence and get the top scoring families (see here: https://pfam-docs.readthedocs.io/en/latest/restful-interface.html)
- Use HMMER webiste to do basically the same, however we can consider other search databases there (https://hmmer-web-docs.readthedocs.io/en/latest/api.html#available-services)
- Run a BLAST and simply parse the output to get the family of the best match

If there are multiple families matching, the output specifies which is the aligned part that matches each family HMM profile.

3. Score each region in the query sequence(s) using the corresponding family profiles, assigning to each AA the corresponding mean normalized b-factor.

## "Training" part:
1. Get the set of reference sequences; should be experimental and from crystallography so that they have b-factors.
**QUESTION: how to get them? we just get N sequences of the corresponding family? **

2. Normalize c-alpha b-factors for each reference protein (rest the mean, divide by std.dev?)

3. Get mean norm-b-factor for each AA in each protein.

4. Get mean of norm-b-factors for each AA, using (weighted?) mean of the intra-protein means. We may use weights based on the quality of the structures (r-factor, resolution: https://proteinstructures.com/structure/protein-databank/) and we may also use only protein structures above a certain quality (although this may be problematic if we get as input some hard-to-crystallyze proteins that have bad structures...)

## First steps to do
These I think we will need to do for sure:
- Define function/class method to assign b-factors to the query
- Define function/method to compute the b-factor profile give a set of model sequences

# DUBTES

* DUBTE 1 -Respecte el càlcul de la flexibilitat quina és la part del aminoàcid que hauríem de tenir en compte (side-chain, backbone, les dues). 

* RESPOSTA 1- Només del backbone, les cadenes laterals es mouen massa I si ho féssiu, hauria de ser  en una conformacio definida.

* DUBTE 2 - Score por aminoacido o por ventana?? 
