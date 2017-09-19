"""Main snakefile"""

include: "const.py"
include: "sf_ann.py"
include: "sf_clinvar.py"
include: "sf_eval.py"
include: "sf_grant.py"
include: "sf_coords.py"

rule all:
    input: DOCS + 'plots/missense_clinvar_roc_feature_union.png',
           DOCS + 'plots/missense_fg_roc.png', WORK + 'eval/missense_fg.dat'

#DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.xls', DATA + 'interim/clinvar/clinvar.dat'
