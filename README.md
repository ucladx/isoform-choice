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
Create a list of ENST for the mane and MSKCC clinical isoforms

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

Wrote a script to find variants where MSK/MANE isoforms disagree on the consequence and amino-acid change.
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

Quantiied levels of actionability per variant for MSKCC and Mane OncoKB annotated documents. This will create two documents, one text file for MANE isoforms and another file for MSKCC isoforms. Can use this to create a bar chart or pie chart in ggplot. 
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
