"""This is a snakefile that runs in docker."""
from const import *

rule plot_gene_missense_counts:
    input:  WORK + 'eval/dat'
    output: DOCS + 'plots/gene_missense_counts.svg'
    shell:  'Rscript {SCRIPTS}plot_gene_counts.R {input} {output}'

rule plot_class_missense_counts:
    input:  WORK + 'eval/dat'
    output: DOCS + 'plots/class_missense_counts.svg'
    shell:  'Rscript {SCRIPTS}plot_benign_path_counts.R {input} {output}'

rule plot_mpc_hist:
    input:  WORK + 'eval/dat'
    output: DOCS + 'plots/mpc_hist.svg'
    shell:  'Rscript {SCRIPTS}plot_mpc_dist.R {input} {output}'

rule all_r_plots:
    input: expand(DOCS + 'plots/{afile}.svg', afile=('mpc_hist', 'class_missense_counts'))
