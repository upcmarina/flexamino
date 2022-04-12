#!/usr/bin/bash

rm all_corr.txt
touch all_corr.txt

for LBL in 1 2 3 4 5 6 7 8 9 10 20 30
do
   out=O08989_$LBL;
   #echo $out
   python3 ../../FlexAmino/flexamino.py -i O08989.fasta -o $out -v -t -r -p $LBL -s 1X1R,1X1S,3KKO,3KKP,3KKQ,3PIR,3PIT;
   
   python3 assess.py -s $LBL;
   cat $out.pearsonr>>all_corr.txt;
done

cat all_corr.txt | sed 's/)/\n/g' | sed 's/(/ /g' | sed 's/,//g'>pearsonR_vs_pbd_lim.txt

python3 corr_vs_PDBlimt.py
