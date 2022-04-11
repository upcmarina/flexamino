**FlexAmino**
==================================

E. Toro, M. Vallejo, S. Vega. **"FlexAmino: assessment of protein flexibility"** Master in Bioinformatics for Health Sciences, 2022

![flexamino_logo](https://user-images.githubusercontent.com/67465839/162766413-d015bd96-3f49-45e0-a2f3-a85dff2070be.png)

# ABOUT

**FlexAmino** is a Python program 
based on the search of homologs for a given protein sequence. Crystallographic beta-factors are extracted from the homologs that have an available known structure and with them, a beta-factor profile is generated and assigned to the query sequence using a structural alignment.

# INSTALLATION

First of all, clone the current repository to your local machine:
```
git clone https://github.com/upcmarina/structural_biology.git
```

Before running **FlexAmino**, you need to install some dependencies. This step can be done by hand or in an automatized way typing the following commands (with root permission):
```
apt-get install python3-setuptools

python3 setup.py install
```

If the installation is successful the following message will appear in the terminal:

![image](https://user-images.githubusercontent.com/67465839/162790974-d571e64e-de0b-4696-b39f-9c2ab6bdcb43.png)

**FlexAmino** uses **clustalw2**. You must copy the bin file provided in the package `bin/clustalw2` to the distribution directory generated in the previous step:
```
cp /bin/clustalw2 /usr/local/lib/python3.8/dist-packages/FlexAmino-1.0-py3.8.egg/EGG-INFO/scripts
```
