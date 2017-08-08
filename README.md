mahdi_epi
==============================

Variant effect prediction using domains

### Methods
* use mutalyzer to convert c. to genome coords (had trouble w/ indels)
* mk vcf file
* annotate with snpeff, dbnsfp
* annotate pfam, hgmd, clinvar, cadd, etc with vxfanno
* ready for gemini db

### variants within pfam domains
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

### pfam domains with more than one variant
* 526 pfam domains hit by at least one variant
* 127 pfam domains with one variant (cannot be evaluated with neighbor variants)

### dist of var counts in pfam domains
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

### pfam domains with the most hits
* 16 Ion_trans:553
* 16 Neur_chan_memb:3
* 18 KCNQ_channel:11
* 18 Neur_chan_memb:129
* 22 IRK:0
* 30 Ion_trans:491
* 53 Neur_chan_memb:71

### var types inside pfam domains
* 4 splice_region_variant
* 54 stop_gained
* 280 synonymous_variant
* 1363 missense_variant

### plan
* fix mutalyzer indels
* find domains with significant burden when compared to exac
* use burden domain feature + other features to predict benign/likely benign vs pathogenic
