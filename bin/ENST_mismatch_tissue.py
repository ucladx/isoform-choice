import argparse, csv
import pandas as pd
import numpy as np
import re

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--msk-enst", action="store", dest="msk", required=True, type=str, help="MSK isoforms with ENSG")
    parser.add_argument("--tissue-enst", action="store", dest="tissue", required=True, type=str, help="Tissue isoforms overrides file")
    parser.add_argument("--tissue-type", action="store", dest="tissuetype", required=True, type=str, help="tissue type of interest (e.g. breast, stomach, lung)")
    args = parser.parse_args()
    
    
    with open(args.msk, 'r') as fh:
        mskdf = pd.read_csv(fh, sep = '\t', header = 0, keep_default_na = False, low_memory = False)
    with open(args.tissue, 'r') as fh:
        tissuedf = pd.read_csv(fh, sep = '\t', header = None, keep_default_na = False, low_memory = False, names = ["ENST_ID", "ENSG_ID","TPM"])
        tissuedf1 = tissuedf.copy()
        tissuedf1['ENSG_ID']=tissuedf['ENSG_ID'].str.replace('\.\d+$','')
        combined= mskdf.merge(tissuedf1, how = "outer", on =('ENSG_ID'), suffixes=('_msk', '_tissue'))
        combined1=(combined[combined['ENST_ID_msk'] != combined['ENST_ID_tissue']])
        combined1.dropna(subset = ["ENST_ID_msk", "ENST_ID_tissue"]).to_csv('data/msk_and_{}_different_ENST.txt'.format(args.tissuetype),sep='\t',header=True, index=False)
if __name__ == "__main__":
    main()
