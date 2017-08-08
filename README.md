mahdi_epi
==============================

Variant effect prediction using domains

### Methods
* use mutalyzer to convert c. to genome coords
* mk vcf file
* annotate with snpeff, dbnsfp, pfam, hgmd, clinvar, cadd, etc

### variants with pfam domain
* 13 BENIGN
* 180 LIKLEY_PATHOGENIC
* 205 LIKELY_BENIGN
* 272 PATHOGENIC
* 1031 VUS

### variants with no pfam domain
* 31 BENIGN
* 88 LIKLEY_PATHOGENIC
* 193 PATHOGENIC
* 394 LIKELY_BENIGN
* 1491 VUS

### domains with more than one variant
* 526 pfam domains hit by at least one variant
* 127 pfam domains with one variant

# dist of var counts in pfam domains
* ex: one domain with 22 variants
* 1 22
* 1 30
* 1 53
* 2 13
* 2 14
* 2 15
* 2 16
* 2 18
* 7 12
* 8 10
* 14 7
* 17 8
* 18 9
* 33 5
* 33 6
* 66 4
* 77 3
* 113 2
