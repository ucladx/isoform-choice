## Study effect of gene isoform choice on clinical variant reporting

Fetch MANE Select transcripts (most prevalent in humans and agreed upon by UCSC/Ensembl):
```
wget -P data https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_0.92/MANE.GRCh38.v0.92.summary.txt.gz
gzip -dc data/MANE.GRCh38.v0.92.summary.txt.gz | awk -F'\t' '($10=="MANE Select") {print}' | cut -f4,8 > data/mane_transcripts.txt
```

Fetch transcripts used by MSKCC for clinical reporting:
```
wget -P data https://raw.githubusercontent.com/mskcc/vcf2maf/bc018db/data/isoform_overrides_at_mskcc
```

::TODO:: Standardize gene names before trying to match them up
