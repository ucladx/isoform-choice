import argparse, csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--msk-maf", action="store", dest="msk", required=True, type=str, help="MSKCC maf annotated from OncoKB")
    parser.add_argument("--mane-maf", action="store", dest="mane", required=True, type=str, help="MANE select maf annotated from OncoKB")
    args = parser.parse_args()

    col_list = ['Hugo_Symbol','Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','HGVSp_Short','Consequence','MUTATION_EFFECT','HIGHEST_LEVEL'] 
    with open(args.mane, 'r') as fh:
        manedf = pd.read_csv(fh, sep = '\t', header = 0, keep_default_na = False, low_memory = False, usecols = col_list)
    with open(args.msk, 'r') as fh1:
        mskdf = pd.read_csv(fh1,sep = '\t', header = 0, keep_default_na = False, low_memory = False, usecols = col_list  )
        combined= mskdf.merge(manedf, how = "outer", on =('Hugo_Symbol','Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2'), suffixes=('_mane', '_msk')) 
        combined[(not combined["HIGHEST_LEVEL_mane"].empty or not combined["HIGHEST_LEVEL_msk"].empty) and combined["HIGHEST_LEVEL_mane"] != combined["HIGHEST_LEVEL_msk"]].to_csv('data/different_highest_level.txt',sep='\t',header=True, index=False)
        
if __name__ == "__main__":
    main()
