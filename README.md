**FlexAmino**
==================================

E. Toro, M. Vallejo, S. Vega. **"FlexAmino: assessment of protein flexibility"** Master in Bioinformatics for Health Sciences, 2022

![flexamino_logo](https://user-images.githubusercontent.com/67465839/162766413-d015bd96-3f49-45e0-a2f3-a85dff2070be.png)

# ABOUT

**FlexAmino** is a Python program based on the search of homologs for a given protein sequence. Crystallographic beta-factors are extracted from the homologs that have an available known structure and with them, a beta-factor profile is generated and assigned to the query sequence using a structural alignment.

# INSTALLATION

In order to install **FlexAmino** first of all you must clone the current repository to your local machine:
```
git clone https://github.com/upcmarina/structural_biology.git
```

Once the repository is cloned, you need to install some dependencies. You have two options, either install them by hand or install them in an automatized way typing the following commands (having root permissions):
```
apt-get install python3-setuptools

python3 setup.py install
```

If the installation is successful the following message will appear in the terminal:

![image](https://user-images.githubusercontent.com/67465839/162790974-d571e64e-de0b-4696-b39f-9c2ab6bdcb43.png)

**FlexAmino** uses clustalw2 via a binary file, which is provided in the repository `bin/clustalw2`. In order to include it in the installation, you must create a directory called `bin` in the distribution package path previously generated, and then copy the binary file:
```
mkdir /usr/local/lib/python3.8/dist-packages/FlexAmino-1.0-py3.8.egg/EGG-INFO/scripts/bin
cp /bin/clustalw2 /usr/local/lib/python3.8/dist-packages/FlexAmino-1.0-py3.8.egg/EGG-INFO/scripts/bin
```

Now you can run `flexamino.py` from any directory in your local machine, remember to have in the directory the fasta file with the fasta sequence you will use as an input:
