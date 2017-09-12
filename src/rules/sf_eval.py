"""Make predictions"""
import pandas

rule clinvar_eval:
    input:  DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.xls',
            DATA + 'interim/clinvar/clinvar.dat'
    output: DOCS + 'plots/missense_clinvar_roc_feature_union.svg',
            WORK + 'eval/auc/missense_clinvar_roc_feature_union',
            WORK + 'eval/missense_clinvar_roc_feature_union.dat'
    shell:  'python {SCRIPTS}clinvar_eval.py {input} {output}'

rule fg_eval:
    input:  DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.xls'
    output: DOCS + 'plots/missense_fg_roc.svg',
            WORK + 'eval/missense_fg.dat'
    shell:  'python {SCRIPTS}eval_fg.py {input} {output}'

rule cat_data:
    input:  WORK + 'eval/missense_fg.dat',
            WORK + 'eval/missense_clinvar_roc_feature_union.dat'
    output: o = WORK + 'eval/dat'
    run:
        df = pandas.concat([pandas.read_csv(f, sep='\t') for f in list(input)])
        df.to_csv(output.o, index=False, sep='\t')

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

rule compose:
    input:  DOCS + 'plots/class_missense_counts.svg',
            DOCS + 'plots/mpc_hist.svg',
            DOCS + 'plots/missense_fg_roc.svg',
            DOCS + 'plots/missense_clinvar_roc_feature_union.svg'
    output: DOCS + 'plots/grant_fig.svg'           
    shell:  'python {SCRIPTS}mk_fig.py {input} {output}'

rule png:
    input:  DOCS + 'plots/grant_fig.svg'
    output: DOCS + 'plots/grant_fig.png'
    shell:  'convert -size 1000x1000 {input} {output}'
