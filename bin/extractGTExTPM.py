'''
@description : Use GTEx TPM file, and Sample Attributes file to find highest median per tissue type
@created :01/05/2021
@author : Eah Alina Keomanee, Cyriac Kandoth
'''

import argparse, gzip, csv, re
from statistics import median


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-attributes", action="store", dest="sample", required=True, type=str, help="GTEx file with sample ID, and tissue")
    parser.add_argument("--GTEx-TPM", action="store", dest="tpm", required=True, type=str, help="Gzipped GTEx file with Sample ID, and TPM data for each transcript")
    parser.add_argument("--tissue-type", action="store", dest="tissueinput", required=True, type=str, help="List the tissue of interest, view the SMTS or SMTSD column on GTEx website for available options (e.g. Breast)")
    args = parser.parse_args()


    tissue_samples = [] #create list of sample_ids for tissue of interest
    with open(args.sample, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            tissue = args.tissueinput
            if tissue in row["SMTSD"] or tissue in row["SMTS"]:
                tissue_samples += [row["SAMPID"].strip()] 
                # print(tissue_samples)
                # exit()
    median_tpm = dict() # List the median TPM per transcript per gene across samples of the specified tissue type 
    with gzip.open(args.tpm, 'rt') as fh:
        reader = csv.DictReader((row for row in fh if row.startswith('transcript_id') or row.startswith('ENST')), delimiter='\t')
        for row in reader:
            gene = row["gene_id"]
            transcript = row["transcript_id"]
            enst_id = re.sub('\.\d+$', '', transcript)
            tpms = [float(row[sample]) for sample in tissue_samples if sample in row]
            if gene not in median_tpm:
                median_tpm[gene]=dict()
            median_tpm[gene][enst_id] = median(tpms)  # e.g  median_tmp[BRAF[ENST0002132145 = 56.8]] , median_tmp[BRAF[ENST00021343756 = 323.8]] 
    for gene in median_tpm.keys():
        most_abundant_transcript = ""
        highest_tpm = 0.0
        for enst_id in median_tpm[gene].keys():
            if highest_tpm < (median_tpm[gene][enst_id]):
                most_abundant_transcript = enst_id
                highest_tpm = median_tpm[gene][enst_id]
        if highest_tpm != 0.0:
            print("\t".join([most_abundant_transcript,gene,str(highest_tpm)]))
if __name__ == "__main__":
    main()
