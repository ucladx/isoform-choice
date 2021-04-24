## Study effect of gene isoform choice on clinical variant reporting

Fetch the single isoform per gene used by Memorial Sloan Kettering Cancer Center (MSKCC) for clinical reporting. To deduplicate CDKN2A's two isoforms, manually select p16INK4a (NM_000077.4) instead of p14ARF (NM_058195.3).
```
curl -sL https://raw.githubusercontent.com/mskcc/vcf2maf/46d276c/data/isoform_overrides_at_mskcc_grch38 | grep -v ^# | grep -wv NM_058195.3 | sort -k2,2 > data/mskcc_clinical_isoforms.txt
```

Fetch MANE Select isoforms (well-supported by experimental data and agreed upon by NCBI/EMBL-EBI). Ignore isoforms tagged "MANE Plus Clinical" for now, since they are subjective choices.
```
curl -sL https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_0.92/MANE.GRCh38.v0.92.summary.txt.gz | gzip -dc | awk -F'\t' 'BEGIN{OFS=FS} ($10=="MANE Select") {print $8,$4,$6}' | sort -k2,2 > data/mane_select_isoforms.txt
```

Find MSK genes that are not in MANE (::TODO:: Find if MANE uses a different gene alias or does not define an isoform for these genes).
```
join -t$'\t' -a1 -j2 -o 1.2,2.2 data/mskcc_clinical_isoforms.txt data/mane_select_isoforms.txt | awk -F'\t' '($2=="") {print}'
```

Find genes where MSKCC and MANE disagree on isoform choice.
```
echo -e "gene_name\tmsk_enst\tmane_enst\tmsk_refseq\tmane_refseq" > data/mskcc_mane_discordant_isoforms.txt
join -t$'\t' -a1 -j2 -o 1.2,1.1,2.1,1.3,2.3 data/mskcc_clinical_isoforms.txt data/mane_select_isoforms.txt | awk -F'\t' 'BEGIN{OFS=FS} ($3!="" && $2!=$3) {print}' >> data/mskcc_mane_discordant_isoforms.txt
```

Downloaded [COSMIC v92 GRCh38 VCF](https://cancer.sanger.ac.uk/cosmic/help/file_download):

```
echo "email@example.com:mycosmicpassword" | base64
```
This will create an authenication string. Example authentication string: ZW1haWxAZXhhbXBsZS5jb206bXljb3NtaWNwYXNzd29yZAo= .  You need to pass the authentication string to the server when you make your request

```
curl -H "Authorization: Basic ZW1haWxAZXhhbXBsZS5jb206bXljb3NtaWNwYXNzd29yZAo=" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v92/VCF/CosmicCodingMuts.vcf.gz

```

That request will return a snippet of JSON containing the link that you need to use to download your file. For example:

```
{
    "url" : "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v92/VCF/CosmicCodingMuts.vcf.gz?AWSAccessKeyId=KFGH85D9KLWKC34GSl88&Expires=1521726406&Signature=Jf834Ck0%8GSkwd87S7xkvqkdfUV8%3D"
}

```
Extract the URL from the JSON response and make another request to that URL to download the file. For example:

```
curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v92/VCF/CosmicCodingMuts.vcf.gz?AWSAccessKeyId=KFGH85D9KLWKC34GSl88&Expires=1521726406&Signature=Jf834Ck0%8GSkwd87S7xkvqkdfUV8%3D" --output cosmic_v92_coding_muts.vcf.gz
```

Subset to [UCLA cancer genes](https://github.com/ucladx/panel-design/blob/master/data/exon_targets_grch38.bed). Ensure that the naming convention between both files are consistent. For example,  UCLA cancer genes first row needs to be updated from chr1 to 1 using the following command:

```
sed 's/^chr//' UCLA_cancer_genes_grch38.bed > UCLA_cancer_genes_nocchr_grch38.bed
```

```
bedtools intersect -a cosmic_v92_coding_muts.vcf.gz -b UCLA_cancer_genes_nochr_grch38.bed -wa > cosmic_v92_UCLA_intersect.vcf
```
Create a list of ENST for the mane and MSKCC clinical isoforms(used vcf2msf version v1.6.19)

```
cut -f1 /src/vcf2maf/data/isoform_overrides_at_mskcc_grch38 | sed '1d;$d' > mskcc_ENST_list_grchr38.txt
cut -f1 /src/vcf2maf/data/isoform_overrides_mane_grch38 | sed '1d;$d' > mane_ENST_list_grchr38.txt
```

Create 2 MAFs, one with variant effects predicted on MSK isoforms and the other on MANE isoforms using the grch38 build:

```
perl ~/src/vcf2maf/vcf2maf.pl  --input-vcf /home/eah19/src/isoform-choice/data/cosmic_v92_UCLA_intersect.vep.vcf --custom-enst data/mskcc_ENST_list_grchr38.txt --ref-fasta /home/eah19/.vep/homo_sapiens/101_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz  --ncbi-build GRCh38 --output-maf /home/eah19/src/isoform-choice/data/cosmic_v92_UCLA_intersect_MSKCC.vep.maf
perl ~/src/vcf2maf/vcf2maf.pl --input-vcf /home/eah19/src/isoform-choice/data/cosmic_v92_UCLA_intersect.vep.vcf --custom-enst data/mane_ENST_list_grchr38.txt --ref-fasta /home/eah19/.vep/homo_sapiens/101_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --ncbi-build GRCh38 --output-maf /home/eah19/src/isoform-choice/data/cosmic_v92_UCLA_intersect_mane.vep.maf
```

Annotate each MAF with clinical actionability from OncoKB:
```
python3 ~/src/oncokb-annotator/MafAnnotator.py -i cosmic_v92_UCLA_intersect_MSKCC.vep.maf -r GRCh38 -b OncoKB token -o cosmic_v92_UCLA_intersect_oncoKB_MSKCC.vep.maf

python3 ~/src/oncokb-annotator/MafAnnotator.py -i cosmic_v92_UCLA_intersect_mane.vep.maf -r GRCh38 -b OncoKB token -o cosmic_v92_UCLA_intersect_oncoKB_mane.vep.maf
```

Find where MSK/MANE isoforms disagree on the consequence and amino-acid changes.
```
python3 bin/compare_isoforms.py --variant-list data/cosmic_v92.ucla_genes.msk_isoforms.oncokb.maf.gz --msk-isoforms data/mskcc_clinical_isoforms.txt --mane-isoforms data/mane_select_isoforms.txt | sort -u > mskcc_mane_discordant_isoforms.txt
```

Created de-deuplicated MSKCC and Mane MAF file. Can also de-duplicate the cosmic_v92_coding_muts.vcf.gz instead of doing this after the annotation steps 

```
    head -n1 cosmic_v92_UCLA_intersect_oncoKB_mane.vep.maf > cosmic_v92_UCLA_intersect_oncoKB_mane.dedup.vep.maf
    tail -n+2 cosmic_v92_UCLA_intersect_oncoKB_mane.vep.maf | sort -u -k5,5V -k6,6n -k11,11 -k13,13 >> cosmic_v92_UCLA_intersect_oncoKB_mane.dedup.vep.maf

    head -n1 cosmic_v92_UCLA_intersect_oncoKB_MSKCC.vep.maf > cosmic_v92_UCLA_intersect_oncoKB_MSKCC.dedup.vep.maf
    tail -n+2 cosmic_v92_UCLA_intersect_oncoKB_MSKCC.vep.maf | sort -u -k5,5V -k6,6n -k11,11 -k13,13 >> cosmic_v92_UCLA_intersect_oncoKB_MSKCC.dedup.vep.maf
```

Quantified levels of actionability per variant for MSKCC and Mane OncoKB annotated documents. This will create two documents, one text file for MANE isoforms and another file for MSKCC isoforms. Can use this to create a bar chart or pie chart in ggplot. 
```
python3 bin/count_conseq.py --msk-maf data/cosmic_v92_UCLA_intersect_oncoKB_MSKCC.dedup.vep.maf --mane-maf data/cosmic_v92_UCLA_intersect_oncoKB_mane.dedup.vep.maf
```
Gene specific levels of actionability for MSKCC and Mane can be quanitified with the following command. Files that are created from this command can be compared. 

```
python3 bin/level_per_gene.py --msk-maf data/cosmic_v92_UCLA_intersect_oncoKB_MSKCC.dedup.vep.maf --mane-maf data/cosmic_v92_UCLA_intersect_oncoKB_mane.dedup.vep.maf 
```
To pinpoint where the differences in levels of actionability for each variant came from, this script prints out where the MSKCC and Mane transcripts do not have the same levels of actionability. 
```
python3 bin/Compare_Highest_Level_Per_Variant.py --msk-maf data/cosmic_v92_UCLA_intersect_oncoKB_MSKCC.dedup.vep.maf --mane-maf data/cosmic_v92_UCLA_intersect_oncoKB_mane.dedup.vep.maf
```
Compare ENST IDs between MSKCC and Mane to create a file where MSKCC and MANE transcripts are discrepant. File ouput is data/msk_and_mane_different_ENST.txt.
```
python3 bin/ENST_mismatch_mane.py --msk-enst data/mskcc_clinical_isoforms.txt --mane-enst data/mane_select_isoforms.txt 

```
Combine variants in one file along with the transcripts, HGVSp notations, and variant consequences from the deduplicated MSKCC and mane vep output maf files. Any empty cells in the HGVSp, variant consequences were dropped into data/non_empty_combined_mutations_mane.txt
```
python3 bin/combine_variant_mane.py --msk-maf data/cosmic_v92_UCLA_intersect_oncoKB_MSKCC.dedup.vep.maf  --tissue-maf data/cosmic_v92_UCLA_intersect_oncoKB_mane.dedup.vep.maf
```
classifying mutation and amino acid differences between MSKCC and MANE
```
python3 bin/classifying_variant_type.py --combined-mutations data/non_empty_combined_mutations_mane.txt > data/variant_type_mane.txt
```

Quantified levels of actionability per variant for MSKCC and Mane OncoKB annotated documents. This will create two documents, one text file for MANE isoforms and another file for MSKCC isoforms. Can use this to create a bar chart or pie chart in ggplot. 
```
python3 bin/count_conseq.py --msk-maf data/cosmic_v92_UCLA_intersect_oncoKB_MSKCC.dedup.vep.maf --mane-maf data/cosmic_v92_UCLA_intersect_oncoKB_mane.dedup.vep.maf
```
Gene specific levels of actionability for MSKCC and Mane can be quanitified with the following command. Files that are created from this command can be compared. 

```
python3 bin/level_per_gene.py --msk-maf data/cosmic_v92_UCLA_intersect_oncoKB_MSKCC.dedup.vep.maf --mane-maf data/cosmic_v92_UCLA_intersect_oncoKB_mane.dedup.vep.maf 
```
To pinpoint where the differences in levels of actionability for each variant came from, this script prints out where the MSKCC and Mane transcripts do not have the same levels of actionability. 
```
python3 bin/Compare_Highest_Level_Per_Variant.py --msk-maf data/cosmic_v92_UCLA_intersect_oncoKB_MSKCC.dedup.vep.maf --mane-maf data/cosmic_v92_UCLA_intersect_oncoKB_mane.dedup.vep.maf
```

AACR Project Genie was used because 1) clinically relevant mutations were identified and 2) cancer type was specificied for the different samples. All mutations in AACR Project Genie were ran in maftovcf
```
perl /home/eah19/bin/vcf2maf/maf2vcf.pl --input-maf data/data_mutations_extended_edited.txt --ref-fasta /home/ckandoth/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz  --per-tn-vcfs  --output-dir /home/eah19/src/isoform-choice/data/genieVCFs

```
Columns 1 through 8 were cut in the VCF file and the file was sorted using bcftools

```
cut -f1-8 data/genieVCFs/data_mutations_extended_edited.vcf | bcftools sort > data/genieVCFs/data_mutations_extended_sorted.vcf

```

The VCF files were run through maf to vcf(v1.6.20) using liftover to convert from a GRch37 to GRCh38 and MSKCC transcripts were used as the custom ENST list. Steps to cut the MSKCC ENST IDs is also listed here
```
cut -f1 /home/eah19/bin/vcf2maf/data/isoform_overrides_at_mskcc_grch38 > /home/eah19/src/isoform-choice/data/mskcc_ENST_list_grchr38

perl /home/eah19/bin/vcf2maf/vcf2maf.pl --input-vcf data/genieVCFs/data_mutations_extended_sorted.vcf --ref-fasta /home/ckandoth/.vep/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --vep-path /home/ckandoth/miniconda3/bin and --vep-data /home/ckandoth/.vep --remap-chain /home/eah19/bin/vcf2maf/data/GRCh37_to_GRCh38.chain --ncbi-build GRCh38 --custom-enst data/mskcc_ENST_list_grchr38 --output-maf /home/eah19/src/isoform-choice/data/genie_merged_grch38_MSKCC.maf 

```

The MSKCC transcript file from AACR Genie Project Data mutations were de-duplicated

```
head -n2 data/genie_merged_grch38_MSKCC.maf > data/genie_merged_grch38_MSKCC_vep_dedup.maf
tail -n+3 data/genie_merged_grch38_MSKCC.maf| sort -u -k5,5V -k6,6n -k11,11 -k13,13 >> data/genie_merged_grch38_MSKCC_vep_dedup.maf

```

GTEX v8 TPM data was used to create lists of the most prevalent transcript in breast, lung and stomach.The transcript with the highest median TPM for each gene was extracted for specificied tissues(can be indicated under --tissue-type). The types of tissues in the GTEx dataset are in the Sample Attributes excel, in the SMTS or SMTSD row. The script will loo to see if tissue-type is specified in the SMTS or SMTSD row and will subset to these samples when calculating the most prevalent transcripts per gene. GTEx data was downloaded frpm here: https://gtexportal.org/home/datasets
```
python3 bin/extractGTExTPM.py --sample-attributes data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt --GTEx-TPM data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz --tissue-type Breast > data/GTEx_Breast_isoform_overrides.txt

```

Next, the data mutations from AACR Project Genie was subsetted to different types of cancers(e.g. Breast,Lung,Sromach). This scripts takes the cancer type and uses it to look through the column cancer type and cancer type detailed. If word is present in these columns, it subsets to these sample IDs. 

```
python3 bin/tissue_specific_genie.py --sample-attributes /home/eah19/src/isoform-choice/data/data_clinical_sample_no_header.txt --data-mutations /home/eah19/src/isoform-choice/data/data_mutations_extended_edited.txt --cancer-type Breast
```
The subsetted file needs to be run through maf to VCF

```
perl /home/eah19/bin/vcf2maf/maf2vcf.pl --input-maf data/data_mutations_extended_breast.txt --ref-fasta /home/ckandoth/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz  --per-tn-vcfs  --output-dir /home/eah19/src/isoform-choice/data/breastgenieVCFs

```

The vcf file is cut from columns 1 through 8 and bcf tools are used to sort the vcf file
```
cut -f1-8 data/breastgenieVCFs/data_mutations_extended_breast.vcf | bcftools sort > data/breastgenieVCFs/data_mutations_extended_breast_sorted.vcf
```
vcf2maf is used, along with liftover to convert GRch37 and CRch38. VEP annotation is also used for the specific transcript IDs that are specified under the custom-ENST list
```
perl /home/eah19/bin/vcf2maf/vcf2maf.pl --input-vcf data/breastgenieVCFs/data_mutations_extended_breast_sorted.vcf --ref-fasta /home/ckandoth/.vep/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --vep-path /home/ckandoth/miniconda3/bin and --vep-data /home/ckandoth/.vep --remap-chain /home/eah19/bin/vcf2maf/data/GRCh37_to_GRCh38.chain --ncbi-build GRCh38 --custom-enst data/GTEx_Breast_isoform_overrides.txt --output-maf /home/eah19/src/isoform-choice/data/genie_merged_grch38_breast.maf 

```

The output maf files were deduplicated
```
head -n2 data/genie_merged_grch38_breast.maf > data/genie_merged_grch38_breast_vep_dedup.maf
tail -n+3 data/genie_merged_grch38_breast.maf | sort -u -k5,5V -k6,6n -k11,11 -k13,13 >> data/genie_merged_grch38_breast_vep_dedup.maf
```
Combine and categorize variant consequences between MSK and Tissue files. Combine variants in one file along with the transcripts, HGVSp notations, and variant consequences from the deduplicated MSKCC and mane vep output maf files. Any empty cells in the HGVSp, variant consequences were dropped into data/non_empty_combined_mutations_mane.txt. Consequences are then categorized into changes in consequence types. 

```
python3 bin/consequence_categorization.py  --msk-maf data/genie_merged_grch38_MSKCC_vep_dedup.maf --tissue-maf data/genie_merged_grch38_breast_vep_dedup.maf --tissue-type breast --msk-genes data/mskcc_gene_names.txt --msk-transcripts data/mskcc_ENST_list_grchr38 --tissue-transcripts data/GTEx_Breast_isoform_overrides.txt > data/consequence_categorization_breast.txt

```
A file to create transcript mismatches between MSKCC transcript and tissue isoform overrides

```
python3 bin/ENST_mismatch_tissue.py --msk-enst data/MSKCC_ENSG.txt --tissue-enst data/GTEx_Breast_isoform_overrides.txt --tissue-type breast

```

The same set of scripts was run for Lung and Stomach
Lung
```
python3 bin/extractGTExTPM.py --sample-attributes data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt --GTEx-TPM data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz --tissue-type Lung > data/GTEx_Lung_isoform_overrides.txt
python3 bin/tissue_specific_genie.py --sample-attributes /home/eah19/src/isoform-choice/data/data_clinical_sample_no_header.txt --data-mutations/home/eah19/src/isoform-choice/data/data_mutations_extended_edited.txt --cancer-type Lung
perl /home/eah19/bin/vcf2maf/maf2vcf.pl --input-maf data/data_mutations_extended_lung.txt --ref-fasta /home/ckandoth/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz  --per-tn-vcfs  --output-dir /home/eah19/src/isoform-choice/data/lunggenieVCFs
cut -f1-8 data/lunggenieVCFs/data_mutations_extended_lung.vcf | bcftools sort > data/lunggenieVCFs/data_mutations_extended_lung_sorted.vcf
perl /home/eah19/bin/vcf2maf/vcf2maf.pl --input-vcf data/lunggenieVCFs/data_mutations_extended_lung_sorted.vcf --ref-fasta /home/ckandoth/.vep/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --vep-path /home/ckandoth/miniconda3/bin and --vep-data /home/ckandoth/.vep --remap-chain /home/eah19/bin/vcf2maf/data/GRCh37_to_GRCh38.chain --ncbi-build GRCh38 --custom-enst data/GTEx_Lung_isoform_overrides.txt --output-maf /home/eah19/src/isoform-choice/data/genie_merged_grch38_lung.maf 
head -n2 data/genie_merged_grch38_lung.maf > data/genie_merged_grch38_lung_vep_dedup.maf
tail -n+3 data/genie_merged_grch38_lung.maf | sort -u -k5,5V -k6,6n -k11,11 -k13,13 >> data/genie_merged_grch38_lung_vep_dedup.maf
python3 bin/consequence_categorization.py  --msk-maf data/genie_merged_grch38_MSKCC_vep_dedup.maf --tissue-maf data/genie_merged_grch38_lung_vep_dedup.maf --tissue-type lung --msk-genes data/mskcc_gene_names.txt --msk-transcripts data/mskcc_ENST_list_grchr38 --tissue-transcripts data/GTEx_Lung_isoform_overrides.txt > data/consequence_categorization_lung.txt
python3 bin/ENST_mismatch_tissue.py --msk-enst data/MSKCC_ENSG.txt --tissue-enst data/GTEx_Lung_isoform_overrides.txt --tissue-type lung

```

Stomach
```
python3 bin/extractGTExTPM.py --sample-attributes data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt --GTEx-TPM data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz --tissue-type Stomach > data/GTEx_Stomach_isoform_overrides.txt
python3 bin/tissue_specific_genie.py --sample-attributes /home/eah19/src/isoform-choice/data/data_clinical_sample_no_header.txt --data-mutations/home/eah19/src/isoform-choice/data/data_mutations_extended_edited.txt --cancer-type Stomach
perl /home/eah19/bin/vcf2maf/maf2vcf.pl --input-maf data/data_mutations_extended_stomach.txt --ref-fasta /home/ckandoth/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz  --per-tn-vcfs  --output-dir /home/eah19/src/isoform-choice/data/stomachgenieVCFs
cut -f1-8 data/stomachgenieVCFs/data_mutations_extended_stomach.vcf | bcftools sort > data/stomachgenieVCFs/data_mutations_extended_stomach_sorted.vcf
perl /home/eah19/bin/vcf2maf/vcf2maf.pl --input-vcf data/stomachgenieVCFs/data_mutations_extended_stomach_sorted.vcf --ref-fasta /home/ckandoth/.vep/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --vep-path /home/ckandoth/miniconda3/bin and --vep-data /home/ckandoth/.vep --remap-chain /home/eah19/bin/vcf2maf/data/GRCh37_to_GRCh38.chain --ncbi-build GRCh38 --custom-enst data/GTEx_Stomach_isoform_overrides.txt --output-maf /home/eah19/src/isoform-choice/data/genie_merged_grch38_stomach.maf
head -n2 data/genie_merged_grch38_stomach.maf > data/genie_merged_grch38_stomach_vep_dedup.maf
tail -n+3 data/genie_merged_grch38_stomach.maf | sort -u -k5,5V -k6,6n -k11,11 -k13,13 >> data/genie_merged_grch38_stomach_vep_dedup.maf
python3 bin/consequence_categorization.py  --msk-maf data/genie_merged_grch38_MSKCC_vep_dedup.maf --tissue-maf data/genie_merged_grch38_stomach_vep_dedup.maf --tissue-type stomach --msk-genes data/mskcc_gene_names.txt --msk-transcripts data/mskcc_ENST_list_grchr38 --tissue-transcripts data/GTEx_Stomach_isoform_overrides.txt > data/consequence_categorization_stomach.txt
python3 bin/ENST_mismatch_tissue.py --msk-enst data/MSKCC_ENSG.txt --tissue-enst data/GTEx_Stomach_isoform_overrides.txt --tissue-type stomach


```
