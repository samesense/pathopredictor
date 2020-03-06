"""Main snakefile"""
include: "const.py"
#include: "sf_patient_counts.py"
#include: "sf_uniprot.py" # including this will grab latest uniprot neutral vars and cause the whole pipeline to rerun
#include: "sf_hgmd.py"

include: "sf_clinvar.py"
include: "sf_ann.py"

FIGS = ('fig1_countPlot',
        'fig3_withinPanel', 'fig4_evalClinvar',
        'fig5_trainClinvarTestPanel',
        'fig6_byGene_and_evalDenovo',)

TABLES = ('S1_missenseDiseaseVariants_hg19', 'S2_missensePredictions_hg19',)

