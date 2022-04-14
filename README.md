**FlexAmino**
==================================

E. Toro, M. Vallejo, S. Vega. **"FLEXAMINO: B-FACTOR PREDICTOR APPROACH FOR THE ASSESSMENT OF PROTEIN FLEXIBILITY"** Master in Bioinformatics for Health Sciences, 2022

![flexamino_logo](https://user-images.githubusercontent.com/67465839/162766413-d015bd96-3f49-45e0-a2f3-a85dff2070be.png)

# ABOUT

**FlexAmino** is a Python program based on the search of homologs for a given protein sequence. Crystallographic beta-factors are extracted from the homologs that have an available known structure and with them, a beta-factor profile is generated and assigned to the query sequence using a structural alignment.

# INSTALLATION

In order to install **FlexAmino** first of all you must clone the current repository to your local machine:
```
git clone https://github.com/upcmarina/flexamino.git
```

Once the repository is cloned, you need to install some dependencies. You have two options, either install them by hand or install them in an automatized way typing the following commands (having root permissions):
```
apt-get install python3-setuptools

python3 setup.py install
```

If the installation is successful the following message will appear in the terminal:

![image](https://user-images.githubusercontent.com/67465839/162790974-d571e64e-de0b-4696-b39f-9c2ab6bdcb43.png)


**FlexAmino** also requires ClustalW to be installed. This can be done using the built-in APT package manager. Then the user only needs to type the following command (having root privileges):
```
apt-get install -y clustalw
```

Now it is possible to run *flexamino.py* from any directory in the local machine. 

# EXECUTION

There are different running options for FlexAmino, depending on the needs of the user.

For first time users, we recommend the **Basic Execution**. With this mode the program is executed from the very beginning to the end, providing the results in a parseable text output and a figure.

But the user may encounter different scenarios where the Basic Execution is not suitable, this is why we provide other options such as recovering BLAST results, among others.

Lastly, please remember that the input files must be in FASTA format with UniProt headers for the program to work properly.

## BASIC EXECUTION

To run the program, simply execute the `flexamino.py` file using the Python3 interpreter and indicate the necessary options. Only the **-i** flag (to indicate the input file) and the **-o** flag (to indicate the prefix of the output files) are mandatory.

For example, we can run the program for the protein sequence Q8IU85. In this case we will also add the **-t** flag (allows us to keep the temporary folder generated during the calculation) and the **-v** flag (for verbose output):
```
flexamino.py -i Q8IU85.fasta -o Q8IU85 -t -v
```

The input file must be in FASTA format with UniProt header. 

The outputname will be used to create the output file names of the text and image files, so no extension should be provided.

## FLAGS

`-h` --help   show this help message and exit

`-i` --input  Mandatory argument. Specify the input FASTA file with the protein sequence to analyze.

`-o` --output   Mandatory argument. Prefix to generate the different output files.

`-v` --verbose  Print the progression of the program execution to the terminal (Standard Error).

`-t` --tmp  Keep the temporary files directory when the program finishes.

`-r` --rescue   Recover a computation from a BLAST result to avoid running BLAST again.

`-p` --pdb_cutoff   Set a maximum number of pdb structures to use for the computation. Default value of 10.

`-w` --winsize   Set a sliding window for smoothing the results. Default value of 1 (i.e. no smoothing)
