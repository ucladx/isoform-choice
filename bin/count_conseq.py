'''
@description : Count the levels of actionability for Onco-KB annotated files
@created : 12/02/2020
@author : Eah Alina Keomanee, Cyriac Kandoth
'''

import argparse, gzip, csv, re
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--msk-maf", action="store", dest="msk", required=True, type=str, help="File with MSKCC isoforms that are annotated with OncoKB")
    parser.add_argument("--mane-maf", action="store", dest="mane", required=True, type=str, help="File with MANE isoforms that are annotated with OncoKB")
    args = parser.parse_args()
    countlevels(args.msk)
    countlevels(args.mane)

def countlevels(filename):
    with open(filename, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t') #goes through header and extracts column names
        level1_count = 0
        level2_count =0
        level3a_count =0
        level3b_count =0
        level4_count = 0
        resistant_count =0
        none_count =0
        for line in reader: #now line is a dictionary aka hashmap that is indexed by column header name
            if line['HIGHEST_LEVEL'] == 'LEVEL_1':
                level1_count +=1 
            elif line['HIGHEST_LEVEL']  == 'LEVEL_2':
                level2_count +=1
            elif line['HIGHEST_LEVEL']  == 'LEVEL_3A':
                level3a_count +=1
            elif line['HIGHEST_LEVEL']  == 'LEVEL_3B':
                level3b_count +=1
            elif line['HIGHEST_LEVEL'] == 'LEVEL_4' :
                level4_count +=1
            elif line['HIGHEST_LEVEL'] == 'LEVEL_R1' or line['HIGHEST_LEVEL'] == 'LEVEL_R2' or line['HIGHEST_LEVEL'] == 'LEVEL_R3':
                resistant_count +=1
            else:
                none_count +=1
    print(level1_count)
    print(level2_count)
    print(level3a_count)
    print(level3b_count)
    print(level4_count)
    print(resistant_count)
    print(none_count)
 
    fields = ['Level', 'Count']
    # data rows of csv file  
    txtfile=filename.replace(".maf",".txt") 
    f = open(txtfile,"w") # opens the csv file with the ability to write into it
    dw = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
    dw.writeheader()
    dw.writerow({'Level': 'Level1', 'Count': level1_count})
    dw.writerow({'Level': 'Level2','Count': level2_count})
    dw.writerow({'Level':'Level3a','Count': level3a_count})
    dw.writerow({'Level':'Level3b','Count': level3b_count})
    dw.writerow({'Level':'Level4', 'Count': level4_count})
    dw.writerow({'Level':'Resistance Count','Count': resistant_count})
    f.close()
if __name__ == "__main__":
    main()