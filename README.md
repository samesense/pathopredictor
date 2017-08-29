mahdi_epi
==============================

Variant effect prediction using domains

### Meeting notes
* [080817](docs/meeting_080817.md)
* [081517](docs/meeting_081517.md)
* [082217](docs/meeting_081517.md)
* [082917](docs/meeting_082917.md)

### Methods
* use mutalyzer to convert c. to genome coords (had trouble w/ indels)
* mk vcf file
* annotate with snpeff, dbnsfp
* annotate 1kg, exac_coverage, evs, kaviar, exac, pfam, hgmd, clinvar, cadd, etc with vcfanno
* limit to rare vars: less than 1% in 1kg total populations
* test for domain enrichment for missense and nonsense
    * gather all rare positions of interest for disease and exac
    * count the number of variant alleles for these positions
    * when a position is not present, assume it is wt for the max number of alleles observed at a position
    * use fisher's exact test for each domain (include no-domain for now)
    * use bh to correct p-values
