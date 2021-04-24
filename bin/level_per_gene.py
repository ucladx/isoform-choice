'''
@description : Count the levels of actionability for each drug for Onco-KB annotated files
@created : 12/04/2020
@author : Eah Alina Keomanee, Cyriac Kandoth
'''

import argparse, gzip, csv, re

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--msk-maf", action="store", dest="msk", required=True, type=str, help="File with MSKCC isoforms that are annotated with OncoKB")
    parser.add_argument("--mane-maf", action="store", dest="mane", required=True, type=str, help="File with MANE isoforms that are annotated with OncoKB")
    args = parser.parse_args()
    msk_counts = countlevels(args.msk)
    mane_counts = countlevels(args.mane)
    print("\t".join(["Gene", "Isoforms", "Level", "Count"]))
    for gene in msk_counts.keys():
        for level in msk_counts[gene].keys():
           print("\t".join([gene, "MSKCC", level, str(msk_counts[gene][level])]))
    for gene in mane_counts.keys():
        for level in mane_counts[gene].keys():
           print("\t".join([gene, "Mane", level, str(mane_counts[gene][level])]))

def countlevels(filename):
    gene_counts = dict()
    with open(filename, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t') #goes through header and extracts column names
        for line in reader:
            gene = line['Hugo_Symbol']
            level = line['HIGHEST_LEVEL'] if line['HIGHEST_LEVEL'] != "" else "None"
            if gene not in gene_counts:
                gene_counts[gene] = dict()
            #gene_counts[gene][level] = gene_counts[gene][level] + 1 if level in gene_counts[gene] else 0
            if level not in gene_counts[gene]:
                gene_counts[gene][level] = 1
            else:
                gene_counts[gene][level] += 1
    return gene_counts

if __name__ == "__main__":
    main()

        