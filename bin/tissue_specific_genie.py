'''
@description : Compare consequences, amino-acid change, and exon number when using different isoforms with Ensembl VEP
@created : 11/11/2020
@author : Eah Alina Keomanee, Cyriac Kandoth
'''

import argparse, csv
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-attributes", action="store", dest="sample", required=True, type=str, help="MSKCC maf annotated from OncoKB")
    parser.add_argument("--data-mutations", action="store", dest="data", required=True, type=str, help="MANE select maf annotated from OncoKB")
    parser.add_argument("--cancer-type", action="store", dest="cancerinput", required=True, type=str, help="List the tissue of interest, view the SMTS column on GTEx website for available options (e.g. Breast)")
    args = parser.parse_args()

    cancer_samples = [] #create list of sample_ids for tissue of interest
    with open(args.sample, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            cancer = args.cancerinput
            if cancer in row["CANCER_TYPE"] or cancer in row["CANCER_TYPE_DETAILED"] :
                cancer_samples += [row["SAMPLE_ID"].strip()] 
    with open(args.data, 'r') as fh:
       mutations = pd.read_csv(fh, sep = '\t', header = 0, keep_default_na = False, low_memory = False)
       mutations[mutations['Tumor_Sample_Barcode'].isin(cancer_samples)].to_csv(('data/data_mutations_extended_stomach.txt'),sep='\t',header=True, index=False)

if __name__ == "__main__":
    main()