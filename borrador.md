**PROJECT INFORMATION**
===================================

Send composition of groups

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

* Use B-factors

* Molecular Dynamics

