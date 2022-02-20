**CALCULATING FLEXIBILITY SCORE**
================================

**MOLT INTERESSANT!!!!!!!!** : Protein flexibility calculation with Python

Use protein flexibility values calculated with the ProtParam module in Biopython (Bio.SeqUtils.ProtParam).  
The flexibility function makes use of amino acid β-factors as reported by Vihinen et al. (1994).  
We were called out in a recent review for (among other things) not using more up-to-date β-factors, particularly those reported in Smith et al. (2003).  
I’m sure there’s some way to wrench open the hood and modify the ProtParam function so that it uses the newer β-factors, but that doesn’t sound like much fun.  
Instead I generated a new function that takes any protein sequence (as a string) and returns a list of flexibilities, smoothed with a 9 residue window.  
I knocked this up quick, so test before you use 

https://www.polarmicrobes.org/protein-flexibility-calculation-with-python/

```
def flexcalc(prot):
    b_factor = {'A': '0.717', 'C': '0.668', 'E': '0.963', 'D': '0.921', 'G': '0.843', 'F': '0.599', 'I': '0.632', 'H': '0.754', 'K': '0.912', 'M': '0.685', 'L': '0.681', 'N': '0.851', 'Q': '0.849', 'P': '0.85', 'S': '0.84', 'R': '0.814', 'T': '0.758', 'W': '0.626', 'V': '0.619', 'Y': '0.615'}
    b_factor_rr = {'A': '0.556', 'C': '0.607', 'E': '0.805', 'D': '0.726', 'G': '0.651', 'F': '0.465', 'I': '0.51', 'H': '0.597', 'K': '0.863', 'M': '0.575', 'L': '0.504', 'N': '0.735', 'Q': '0.729', 'P': '0.752', 'S': '0.698', 'R': '0.676', 'T': '0.648', 'W': '0.577', 'V': '0.503', 'Y': '0.461'}
    b_factor_rf = {'A': '0.704', 'C': '0.671', 'E': '0.911', 'D': '0.889', 'G': '0.811', 'F': '0.582', 'I': '0.617', 'H': '0.734', 'K': '0.862', 'M': '0.641', 'L': '0.65', 'N': '0.848', 'Q': '0.817', 'P': '0.866', 'S': '0.846', 'R': '0.807', 'T': '0.742', 'W': '0.609', 'V': '0.603', 'Y': '0.567'}
    b_factor_ff = {'A': '0.847', 'C': '0.67', 'E': '1.11', 'D': '1.055', 'G': '0.967', 'F': '0.653', 'I': '0.686', 'H': '0.894', 'K': '1.016', 'M': '0.74', 'L': '0.788', 'N': '0.901', 'Q': '1.007', 'P': '0.857', 'S': '0.914', 'R': '0.942', 'T': '0.862', 'W': '0.656', 'V': '0.707', 'Y': '0.741'}
    
    rigid = set(['A', 'C', 'F', 'I', 'H', 'M', 'L', 'W', 'V', 'Y'])
    flex = set(['E', 'D', 'G', 'K', 'N', 'Q', 'P', 'S', 'R', 'T'])
    
    ## get beta factors for each residue, according to neighbors
    
    b_factors = []
                    
    for r in range(0, len(prot)):
        b_query = False ## just making sure something doesn't go awry and an old value is used
        query = prot[r]
        if r == 0:
            b_query = b_factor[query]
        elif r == len(prot) - 1:
            b_query = b_factor[query]
        else:
            if prot[r - 1] in rigid:
                if prot[r + 1] in rigid:
                    b_query = b_factor_rr[query]
                elif prot[r + 1] in flex:
                    b_query = b_factor_rf[query]
            elif prot[r - 1] in flex:
                if prot[r + 1] in rigid:
                    b_query = b_factor_rf[query]
                elif prot[r + 1] in flex:
                    b_query = b_factor_ff[query]
        if b_query == False:
            print(r, 'bad b-factor!')
            exit()
        else:
            b_factors.append(b_query)
            
    ## average over 9 aa window
            
    b_factors_win = []
            
    for b in range(0, len(b_factors)):
        win = b_factors[b:b + 9]
        if len(win) == 9:
            b_win = sum(map(float, win)) / 9.0
            b_factors_win.append(b_win)
            
    return(b_factors_win)
```

Test:

```
prot = 'MKHFSKLCFLLSTFAVSIAPVTWAHEGATHQHANVSKLTDAYTYANYDQVKATHVYLDLNVDFDKKSLSG'\
'FAELSLDWFTDNKAPLILDTRDLVIHRVMAKNSQGQWVKVNYDLAKRDDVLGSKLTINTPLNAKKVRVYY'\
'NSTEKATGLQWLSAEQTAGKEKPFLFSQNQAIHARSWIPIQDTPSVRVTYTARITTDKDLLAVMSANNEP'\
'GTERDGDYFFSMPQAIPPYLIAIGVGDLEFKAMSHQTGIYAESYILDAAVAEFDDTQAMIDKAEQMYGKY'\
'RWGRYDLLMLPPSFPFGGMENPRLSFITPTVVAGDKSLVNLIAHELAHSWSGNLVTNESWRDLWLNEGFT'\
'SYVENRIMEAVFGTDRAVMEQALGAQDLNAEILELDASDTQLYIDLKGRDPDDAFSGVPYVKGQLFLMYL'\
'EEKFGRERFDAFVLEYFDSHAFQSLGTDNFVKYLKANLTDKYPNIVSDNEINEWIFKAGLPSYAPQPTSN'\
'AFKVIDKQINQLVTDELTLEQLPTAQWTLHEWLHFINNLPVDLDHQRMVNLDKAFDLTNSSNAEIAHAWY'\
'LLSVRADYKEVYPAMAKYLKSIGRRKLIVPLYKELAKNAESKAWAVEVYKQARPGYHGLAQGTVDGVLK'

test1 = flexcalc(prot)

from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
test_seq = PA(prot)
test2 = PA.flexibility(test_seq)

import matplotlib.pyplot as plt

plt.scatter(test1[1:], test2)
```

ProtParam values on the y-axis and our values on the x-axis.  Note that the scales are different, they are arbitrary.  I haven’t figured out how to fit a linear model in Python yet (on the to do list), but the fit looks pretty good.  Plenty of scatter but the methods are in general agreement (which is good – we already have one paper published using the old method!).  Time to go back and redo a lot of work from the last two years…

![image](https://user-images.githubusercontent.com/67465839/154842241-26105243-8a67-47e9-8519-a1f32963d50d.png)
