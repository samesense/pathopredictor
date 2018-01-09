"""Predict stats for vars in clinvar and denovo-db"""
import pandas as pd
include: "const.py"

from snakemake.utils import R

rule eval_clinvar_global:
    input:  DATA + 'interim/{dat}/{dat}.limit3.dat',
            DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/uc.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/other/other.eff.dbnsfp.anno.hHack.dat.limit.xls'
    output: WORK + 'global.eval_{dat}.stats',
            WORK + 'global.eval_{dat}.eval'
    shell:  'python {SCRIPTS}score_other_global_model.py {input} {output}'

rule limit_for_plot:
    input:  WORK + '{method}.eval_{dat}.eval'
    output: WORK + '{method}.eval_{dat}.totWrong'
    shell:  "grep 'TotWrong\|count' {input} | grep -v ssue | grep -v earing > {output}"

rule cat:
    input:  expand( WORK + '{{method}}.eval_{dat}.totWrong', dat=('clinvar', 'denovo', 'clinvar_mult', 'clinvar_single', 'clinvar_exp') )
    output: o=WORK + 'totWrong/{method}'
    run:
        pd.concat( [pd.read_csv(x, sep='\t') for x in list(input)] ).to_csv(output.o, index=False, sep='\t')
        
rule plot:
    input:  WORK + 'totWrong/{method}'
    output: DOCS + 'plot/{method}.other.totWrong.png'
    run:
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          p = ggplot(data=d) +
          geom_col(aes(y=var_count,x=score_type, fill=score_type)) +
          facet_grid(clinvar_type~., scale='free') + theme_bw() +
          ylab('Wrong Predictions') +
          theme(axis.text.x = element_text(angle=45, hjust=1)) +
          xlab('') + theme(legend.position="none")
          ggsave("{output}", p)
          """)

rule all_eval:
    input: expand( DOCS + 'plot/{method}.other.totWrong.png', method=('global',) )
    
