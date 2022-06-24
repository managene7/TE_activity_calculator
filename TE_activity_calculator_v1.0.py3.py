# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 15:29:16 2021

@author: minkj
"""

#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-min_score':60,'-fasta':""}
for i in range(len(args)):
    if args[i].startswith("-"):
        try:
            option_dict[args[i]]=args[i+1]
        except:
            if args[0]=="-help":
                print("""
_________________________________________________________________________________

Usage;

-fasta          file name of TE paralogs (DNA sequence)
-min_score      minimum score for blast result parsing (option, default: 60)

_________________________________________________________________________________
""")
                quit()
                
if option_dict["-fasta"]=="":
    print("""
    ___________________________________________________________________
    
    Please enter file name.
    
    To see options;
    python TE_activity_calculator.py -help
    
    Example;
    python TE_activity_calculator.py -fasta test_sequence.txt -min_score 60
    
    ____________________________________________________________________
    
    
    """)
    quit()

import os

current_directory=os.getcwd()

print("Running blast search..")
blast=os.system("blastn -task blastn -outfmt 5 -query %s -subject %s -out %s -max_target_seqs 5000" % (option_dict['-fasta'],option_dict['-fasta'],option_dict['-fasta']+".xml" ))
if blast==0:
    parsing=os.system("python3 %s\\codes\\blast_parsing___xml_format_v2.6.py3.py -ov_len 10000000 -lim_score %s -num_hsp 1 -lim_evalue 1 -xml %s" % (current_directory, option_dict['-min_score'], option_dict['-fasta']+".xml" ))
    if parsing==0:
        activity=os.system("python3 %s\\codes\\TE_activity_calculation_by_UPGMA_v.1.0.py3.py -csv %s" % (current_directory, option_dict['-fasta']+".csv" ))

print("Completed!")
