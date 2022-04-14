**FlexAmino**
==================================

E. Toro, M. Vallejo, S. Vega. **"FLEXAMINO: B-FACTOR PREDICTOR APPROACH FOR THE ASSESSMENT OF PROTEIN FLEXIBILITY"** Master in Bioinformatics for Health Sciences, 2022

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


**FlexAmino** also requires ClustalW to be installed. This can be done using the built-in APT package manager. Then the user only needs to type the following command (having root privileges):
```
apt-get install -y clustalw
```

Now it is possible to run *flexamino.py* from any directory in the local machine. 
