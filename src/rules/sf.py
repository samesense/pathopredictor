"""Main snakefile"""
# import pandas as pd
# from functools import reduce
# from itertools import combinations, chain
# from sklearn import metrics
# from snakemake.utils import R
# from collections import defaultdict
# from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
# from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
# FTP = FTPRemoteProvider()
# from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
# HTTP = HTTPRemoteProvider()
# import os, sys
# from p_change import *

include: "const.py"
include: "sf_patient_counts.py"
#include: "sf_uniprot.py"
include: "sf_hgmd.py"

include: "sf_clinvar.py"
include: "sf_ann.py"
include: "sf_other_disease.py"
include: "sf_gnomad.py"
include: "sf_ndenovo.py"

include: "counting_funcs.py"
include: "calc_wrong_funcs.py"
include: "sf_eval_panel.py"
include: "improve_prob.py"
include: "sf_coords.py"
include: "sf_eval_ahmad.py" # plot_ahmad
include: "sf_gene_pr_curve.py"
include: "sf_paper_data.py"
include: "paper.plots.py"
include: "sf_predict_clinvar.py"
include: "sf_feature_importances.py"
include: "sf_feature_cor.py"

include: "sf.rank.eval.py"
include: "sf_eval_roc.py"

#include: "sf.webtool.py"
#include: "sf_single_gene.py"

rule all_dat:
    input: DATA + 'interim/denovo/denovo.limit3.dat', \
           DATA + 'interim/clinvar/clinvar.limit3.dat', \
           DATA +  'interim/other/other.eff.dbnsfp.anno.hHack.dat.limit.xls',  \
           expand(DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.hHack.dat.limit.xls', lab=('uc', 'EPIv6') ), \
           TMP + 'trees/revel-ccr-is_domain', \
           TMP + 'trees/revel-is_domain', \
           TMP + 'trees/revel', \

FIGS = ('fig1_countPlot', 'fig2_featureImportance', 'fig3_featureCor',
        'fig5_panelEval', 'fig7_byGene_and_evalDenovo', 'fig6_evalClinvar')

TABLES = ('S1_trainingData_hg19',)

rule all_paper_plots:
    input: expand(DOCS + 'paper_plts/{fig}.pdf', fig=FIGS)

rule upload_all:
    input: expand(DBox.remote('ahmad_predictor/{fig}.tiff'), fig=FIGS),
           expand(DBox.remote('ahmad_predictor/{table}.csv'), table=TABLES)
