'''
@description : Compare consequences, amino-acid change, and exon number when using different isoforms with Ensembl VEP
@created : 11/11/2020
@author : Eah Alina Keomanee, Cyriac Kandoth
'''

import argparse, gzip, csv, re

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant-list", action="store", dest="vars", required=True, type=str, help="bgzipped MAF-like list of variants with all_effects column from vcf2maf")
    parser.add_argument("--msk-isoforms", action="store", dest="msk", required=True, type=str, help="list of transcript ENST IDs from MSKCC")
    parser.add_argument("--mane-isoforms", action="store", dest="mane", required=True, type=str, help="list of transcript ENST IDs from MANE Select")
    args = parser.parse_args()

    # Load the MSK/MANE isoform ENST IDs into dictionaries, removing version numbers in the process
    msk_isoform = dict()
    mane_isoform = dict()
    with open(args.msk, 'r') as fh:
        for line in fh:
            enst_id, gene, refseq_id = line.rstrip().split("\t")
            enst_id = re.sub('\.\d+$', '', enst_id)
            msk_isoform[enst_id] = gene
    with open(args.mane, 'r') as fh:
        for line in fh:
            enst_id, gene, refseq_id = line.rstrip().split("\t")
            enst_id = re.sub('\.\d+$', '', enst_id)
            mane_isoform[enst_id] = gene

    # For each variant, compare "Consequence" between MSK/MANE isoforms
    with gzip.open(args.vars, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for variant in reader:
            var_gene = variant['Hugo_Symbol']
            msk_consequence = ""
            mane_consequence = ""
            for effect in variant['all_effects'].split(";"):
                # Skip empty effects caused by the ";" that vcf2maf inserts at the end of the all_effects values
                if(effect != ""):
                    effect_gene, consequence, aa_change, enst_id, refseq_ids = effect.split(",", 4)
                    if(enst_id in msk_isoform and msk_isoform[enst_id] == var_gene):
                        msk_consequence = ",".join([effect_gene, consequence, aa_change])
                    if(enst_id in mane_isoform and mane_isoform[enst_id] == var_gene):
                        mane_consequence = ",".join([effect_gene, consequence, aa_change])
            if(msk_consequence != "" and mane_consequence != "" and msk_consequence != mane_consequence):
                print("\t".join([msk_consequence, mane_consequence]))

if __name__ == "__main__":
    main()
