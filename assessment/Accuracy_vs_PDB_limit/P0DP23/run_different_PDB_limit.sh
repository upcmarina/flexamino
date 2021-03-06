#!/usr/bin/bash

rm all_corr.txt
touch all_corr.txt

for LBL in 1 2 3 4 5 6 7 8 9 10 20 30
do
   out=P0DP23_$LBL;
   #echo $out
   python3 ../../../FlexAmino/flexamino.py -i P0DP23.fasta -o $out -v -t -r -p $LBL -s 1CDL,1CLL,1CTR,1IWQ,1K90,1K93,1L7Z,1LVC,1PK0,1S26,1SK6,1UP5,1WRZ,1XFU,1XFV,1XFW,1XFX,1XFY,1XFZ,1Y6W,1YR5,1YRT,1YRU,1ZOT,1ZUZ,2BE6,2F3Y,2F3Z,2I08,2R28,2V01,2V02,2VAY,2W73,2WEL,2X0G,2Y4V,3BYA,3DVE,3DVJ,3DVK,3DVM,3EWT,3EWV,3G43,3HR4,3O77,3O78,3OXQ,3SUI,3UCT,3UCW,3UCY,4BW7,4BW8,4BYF,4DCK,4DJC,4GOW,4JPZ,4JQ0,4L79,4LZX,4M1L,4OVN,4Q57,4Q5U,4UMO,4UPU,4V0C,5COC,5DBR,5DOW,5DSU,5I0I,5J03,5JQA,5JTH,5K8Q,5V02,5V03,5V7X,5WBX,5WC5,6B8L,6B8M,6B8N,6B8P,6B8Q,6DAD,6DAE,6DAF,6DAH,6EEB,6HCS,6HR1,6K4K,6K4L,6K4R,6M7H,6MUD,6MUE,6N5W,6O5G,6OS4,6PAW,6U39,6U3A,6U3B,6U3D,6XXX,6XY3,6XYR,6Y4P,6YNS,6YNU,7BF1,7BF2,7KL5;
   
   python3 assess.py -s $LBL;
   cat $out.pearsonr>>all_corr.txt;
done

cat all_corr.txt | sed 's/)/\n/g' | sed 's/(/ /g' | sed 's/,//g'>pearsonR_vs_pbd_lim.txt

python3 corr_vs_PDBlimt.py
