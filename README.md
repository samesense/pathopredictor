PathoPredictor
==============================

Missense variant effect pathogenicity prediction using disease specific classifiers for cardiomyopathy, epilepsy, and rasopathies. If you are here to get PathoPredictor scores for a vcf file, see [prediction instructions](docs/how_to_predict.md). If you want to recreate the whole study, see [reproducing](docs/reproducing.md).

### Methods
* use mutalyzer to convert c. to genome coords (had trouble w/ indels)
* mk vcf file
* annotate with snpeff, dbnsfp
* annotate 1kg, exac_coverage, evs, kaviar, exac, pfam, hgmd, clinvar, cadd, etc with vcfanno
