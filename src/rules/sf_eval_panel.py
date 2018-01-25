"""Predict status for gene panel vars"""
include: "const.py"

from snakemake.utils import R

rule eval_panel_global:
    input:  expand(DATA + 'interim/{dat}/{dat}.limit3.dat', dat=('clinvar', 'clinvar_single', 'clinvar_mult', 'clinvar_exp', 'denovo')),
            DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/uc.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/other/other.eff.dbnsfp.anno.hHack.dat.limit.xls'
    output: WORK + 'global.eval_panel.{cols}.stats',
            WORK + 'global.eval_panel.{cols}.eval'
    shell:  'python {SCRIPTS}score_panel_global_model.py {wildcards.cols} {input} {output}'

rule eval_panel_single_gene:
    input:  expand(DATA + 'interim/clinvar{dat}/{dat}.limit3.dat', dat=('clinvar', 'clinvar_single', 'clinvar_mult', 'clinvar_exp', 'denovo')),
            DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/uc.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/other/other.eff.dbnsfp.anno.hHack.dat.limit.xls'
    output: WORK + 'single.eval_panel.stats',
            WORK + 'single.eval_panel.eval'
    shell:  'python {SCRIPTS}score_panel_single_gene_model.py {input} {output}'
    
rule limit_for_plot:
    input:  WORK + '{method}.eval_panel.{cols}.eval'
    output: WORK + '{method}.eval_panel.{cols}.totWrong'
    shell:  "grep 'TotWrong\|count' {input} | grep -v ssue | grep -v earing > {output}"

rule plot:
    input:  WORK + '{method}.eval_panel.{cols}.totWrong'
    output: DOCS + 'plot/{method}.eval_panel.{cols}.totWrong.png'
    run:
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          p = ggplot(data=d) +
          geom_col(aes(y=var_count,x=score_type, fill=score_type)) +
          facet_grid(disease~., scale='free') + theme_bw() +
          ylab('Wrong Predictions') +
          theme(axis.text.x = element_text(angle=45, hjust=1)) +
          xlab('') + theme(legend.position="none")
          ggsave("{output}", p)
          """)

rule all_eval:
    input: expand( DOCS + 'plot/{method}.eval_panel.{cols}.totWrong.png', method=('global',), cols=('mpc', 'revel', 'mpc-revel', 'ccr', 'mpc-revel-ccr', 'mpc-ccr', 'revel-mpc') )

# ggplot(data=d) + geom_col(aes(y=var_count,x=score_type, fill=score_type)) + facet_grid(disease~., scale='free') + theme_bw() + ylab('Wrong Predictions') + theme(axis.text.x = element_text(angle=90)) + xlab('') + theme(legend.position="none")    


