import argparse, csv
import numpy as np
import pandas as pd

#maybe I can use subprocess for the command line vep/bcftools steps

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--msk-maf", action="store", dest="msk", required=True, type=str, help="maf created using preferred transcript list that is annotated with vep")
    parser.add_argument("--tissue-maf", action="store", dest="tissue", required=True, type=str, help="maf created using transcripts using most prevalent transcripts in different tissues that is annotated with vep")
    parser.add_argument("--tissue-type", action="store", dest="tissuetype", required=True, type=str, help="tissue type of interest (e.g. breast, stomach, lung), will be used to name the file that is combined with preferred transcripts")
    parser.add_argument("--msk-genes", action="store", dest="genes", required=True, type=str, help="List of gene names used in this analysis")
    parser.add_argument("--msk-transcripts", action="store", dest="msktransc", required=True, type=str, help="Preferred transcript IDs that will be used to double check the correct preferred transcript ID")
    parser.add_argument("--tissue-transcripts", action="store", dest="tissuetransc", required=True, type=str, help="Tissue transcript IDs that will be used to double check the correct transcript ID")
    args = parser.parse_args()
    print("\t".join(["Location","Hugo_Symbol","Variant_Classification_tissue","Transcript_ID_tissue","HGVSp_tissue","Variant_Classification_msk","Transcript_ID_msk","HGVSp_msk","Consequence Category"]))

    col_list = ['Hugo_Symbol','Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Transcript_ID','HGVSp_Short','Variant_Classification'] 
    with open(args.tissue, 'r') as fh:
        tissuemaf = pd.read_csv(fh, sep = '\t', header = 1, keep_default_na = False, low_memory = False, usecols = col_list)
    #consequence_type = []
    with open(args.msk, 'r') as fh:
        mskmaf = pd.read_csv(fh,sep = '\t', header = 1, keep_default_na = False, low_memory = False, usecols = col_list )
        combined= mskmaf.merge(tissuemaf, how = "outer", on =('Hugo_Symbol','Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2'), suffixes=('_msk', '_tissue')) 
        #df.to_csv('myfile_{}.csv'.format(gene_name), sep=',', index=False)
        combined.dropna(subset = ["Variant_Classification_tissue", "Variant_Classification_msk"]).to_csv('data/non_empty_combined_mutations_{}.txt'.format(args.tissuetype),sep='\t',header=True, index=False)

    variant_classification = dict()
    with open('data/non_empty_combined_mutations_{}.txt'.format(args.tissuetype),'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            chromosome = row["Chromosome"]
            start = row["Start_Position"]
            ref_allele = row["Reference_Allele"]
            alt_allele = row["Tumor_Seq_Allele2"]
            unique = chromosome+":"+start+":"+ref_allele+":"+alt_allele
            HGVSp_tissue= str(row["HGVSp_Short_tissue"].strip())
            HGVSp_msk = str(row["HGVSp_Short_msk"].strip())
            variant_classification_tissue = str(row["Variant_Classification_tissue"].strip())
            variant_classification_msk = str(row["Variant_Classification_msk"].strip())
            if unique not in variant_classification:
                variant_classification[unique] = dict()
                truncating =  ["Frame_Shift_Del" , "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation"]
                non_truncating = ["3'Flank","3'UTR", "5'Flank", "5'UTR","Intron","Missense_Mutation", "Nonstop_Mutation","RNA", "Silent", "In_Frame_Del", "In_Frame_Ins","IGR","Splice_Region"]
                non_truncating_coding = ["Missense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins"]
                non_truncating_non_coding = ["3'Flank","3'UTR", "5'Flank", "5'UTR","Intron","IGR"]
                non_truncating_silent = ["Silent"]
                non_truncating_rna_splice_region = ["RNA", "Splice_Region"]
                if row["Variant_Classification_tissue"] in truncating and row["Variant_Classification_msk"] in non_truncating: 
                    variant_classification[unique] = "MP:nontruncating to truncating"
                elif row["Variant_Classification_tissue"] in non_truncating and row["Variant_Classification_msk"] in truncating: 
                    variant_classification[unique] = "MB:truncating to nontruncating"
                elif HGVSp_tissue == HGVSp_msk and variant_classification_tissue == variant_classification_msk:
                    variant_classification[unique] = "same"
                elif HGVSp_tissue != HGVSp_msk and variant_classification_tissue == variant_classification_msk:
                    variant_classification[unique] = "Likely No Change"
                elif variant_classification_tissue in non_truncating_coding and variant_classification_msk in non_truncating_non_coding:
                    variant_classification[unique] = "MP: Non Coding to Coding" 
                elif variant_classification_tissue in non_truncating_non_coding and variant_classification_msk in non_truncating_coding:
                    variant_classification[unique] = "MB: Coding to Non Coding" 
                elif variant_classification_tissue in non_truncating_silent  and  variant_classification_msk in non_truncating_coding:
                    variant_classification[unique] = "MB: Coding to Silent" 
                elif variant_classification_tissue in non_truncating_coding and variant_classification_msk in non_truncating_silent:
                    variant_classification[unique] = "MP: Silent to Coding" 
                elif  variant_classification_tissue in non_truncating_rna_splice_region  and variant_classification_msk in non_truncating_non_coding:
                    variant_classification[unique] = "MP: Non Coding to RNA/Splice region" 
                elif variant_classification_tissue in non_truncating_non_coding  and variant_classification_msk in non_truncating_rna_splice_region:
                    variant_classification[unique] = "MB: RNA/Splice Region to Non Coding" 
                elif  variant_classification_tissue in non_truncating_rna_splice_region  and variant_classification_msk in non_truncating_silent:
                    variant_classification[unique] = "MP: Silent to RNA/Splice region" 
                elif variant_classification_tissue in non_truncating_silent and variant_classification_msk in non_truncating_rna_splice_region:
                    variant_classification[unique] = "MB: RNA/Splice Region to Silent" 
                elif (variant_classification_tissue in non_truncating_non_coding or variant_classification_tissue in non_truncating_rna_splice_region) and (variant_classification_msk in non_truncating_non_coding or variant_classification_msk in non_truncating_rna_splice_region):
                    variant_classification[unique] = "Likely No Change"
                elif variant_classification_tissue in truncating and variant_classification_msk in truncating:
                    variant_classification[unique] = "Likely No Change"
                elif (variant_classification_tissue in non_truncating_non_coding or variant_classification_tissue in non_truncating_silent) and (variant_classification_msk in non_truncating_non_coding or variant_classification_msk in non_truncating_silent):
                    variant_classification[unique] = "Likely No Change"
                elif (variant_classification_tissue in non_truncating_rna_splice_region or variant_classification_tissue in non_truncating_coding) and (variant_classification_msk in non_truncating_rna_splice_region or variant_classification_msk in non_truncating_coding):
                    variant_classification[unique] = "Likely No Change"
                else: 
                    variant_classification[unique] = "everything else"
                classifications = ["MP:nontruncating to truncating", "MB:truncating to nontruncating","Likely No Change","MP: Non Coding to Coding","MB: Coding to Non Coding" , "MB: Coding to Silent" , "MP: Silent to Coding","MP: Non Coding to RNA/Splice region", "MB: RNA/Splice Region to Non Coding","everything else"]
                genes = []
                MSK_transcript = []
                Tissue_transcript = []
                with open(args.genes, 'r') as fh:
                    for line in fh:
                        genes += [line.strip()]
                with open(args.msktransc, 'r') as fh:
                    for line in fh:
                        MSK_transcript += [line.strip()]
                        # print(MSK_transcript)
                with open(args.tissuetransc, 'r') as fh:
                    Tissue_transcript = [line.split()[0] for line in fh]
                    # print(Tissue_transcript)
                    if (variant_classification[unique]in classifications) and row["Hugo_Symbol"] in genes and row["Transcript_ID_tissue"] in Tissue_transcript and row["Transcript_ID_msk"] in MSK_transcript :
                        if row["Transcript_ID_msk"] not in ['ENST00000579755', 'ENST00000304494'] and row["Transcript_ID_tissue"] not in ['ENST00000579755', 'ENST00000304494']:
                            print("\t".join([unique,row["Hugo_Symbol"], row["Variant_Classification_tissue"], row["Transcript_ID_tissue"],HGVSp_tissue, row["Variant_Classification_msk"], row["Transcript_ID_msk"], HGVSp_msk, variant_classification[unique]]))
                   
if __name__ == "__main__":
    main()
