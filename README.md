## Study effect of gene isoform choice on clinical variant reporting

Fetch the single isoform per gene used by Memorial Sloan Kettering Cancer Center (MSKCC) for clinical reporting. To deduplicate CDKN2A's two isoforms, manually select p16INK4a (NM_000077.4) instead of p14ARF (NM_058195.3).
```
curl -sL https://raw.githubusercontent.com/mskcc/vcf2maf/46d276c/data/isoform_overrides_at_mskcc_grch38 | grep -v ^# | grep -wv NM_058195.3 | sort -k2,2 > data/mskcc_clinical_isoforms.txt
```

Fetch MANE Select isoforms (well-supported by experimental data and agreed upon by NCBI/EMBL-EBI). Ignore isoforms tagged "MANE Plus Clinical" for now, since they are subjective choices.
```
curl -sL https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_0.92/MANE.GRCh38.v0.92.summary.txt.gz | gzip -dc | awk -F'\t' 'BEGIN{OFS=FS} ($10=="MANE Select") {print $8,$4,$6}' | sort -k2,2 > data/mane_select_isoforms.txt
```

::TODO:: Standardize gene names before trying to match them up

Find genes where MSKCC and MANE disagree on isoform choice.
```
echo -e "gene_name\tmsk_enst\tmane_enst\tmsk_refseq\tmane_refseq" > data/mskcc_mane_discordant_isoforms.txt
join -t$'\t' -a1 -j2 -o 1.2,1.1,2.1,1.3,2.3 data/mskcc_clinical_isoforms.txt data/mane_select_isoforms.txt | awk -F'\t' 'BEGIN{OFS=FS} ($3!="" && $2!=$3) {print}' >> data/mskcc_mane_discordant_isoforms.txt
```
