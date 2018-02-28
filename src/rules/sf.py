"""Main snakefile"""
import pandas as pd
from functools import reduce
from itertools import combinations, chain
from snakemake.utils import R

include: "const.py"
include: "sf_ann.py"
include: "sf_other_disease.py"
include: "sf_clinvar.py"
include: "sf_denovo.py"

include: "sf_eval_panel.py"
include: "sf_coords.py"
include: "sf_eval_ahmad.py" # plot_ahmad
include: "paper.plots.py"
include: "sf_predict_clinvar.py"

include: "sf.webtool.py"

rule all_dat:
    input: DATA + 'interim/denovo/denovo.limit3.dat', \
           DATA + 'interim/clinvar/clinvar.limit3.dat', \
           DATA +  'interim/other/other.eff.dbnsfp.anno.hHack.dat.limit.xls',  \
           expand(DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.hHack.dat.limit.xls', lab=('uc', 'EPIv6') ), \
           TMP + 'trees/revel-ccr-is_domain', \
           TMP + 'trees/revel-is_domain', \
           TMP + 'trees/revel', \





