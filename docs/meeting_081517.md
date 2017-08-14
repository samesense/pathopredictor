## 081517

### Methods
* use mutalyzer to convert c. to genome coords (had trouble w/ indels, the hg19 coords in xlsx are wrong)
* mk vcf file
* annotate with snpeff, dbnsfp
* annotate 1kg, exac_coverage, evs, kaviar, exac, pfam, hgmd, clinvar, cadd, etc with vxfanno
* ready for gemini db
* have exac in gemini db

### fg summary

### exac summary

### variant positions within pfam domains
* 13 BENIGN
* 180 LIKLEY_PATHOGENIC
* 205 LIKELY_BENIGN
* 272 PATHOGENIC
* 1031 VUS

### pfam domains with the most hits (all clinical classifications)
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
* get read genome coords?
* fix mutalyzer indels
* compare synon vars to exac to find correction factor
* limit to rare vars? and find domains with significant burden when compared to exac
* use burden domain feature + other features to predict benign/likely benign vs pathogenic
