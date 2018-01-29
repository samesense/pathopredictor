"""Predict status for gene panel vars"""
import pandas as pd
from functools import reduce
from itertools import combinations, chain

include: "const.py"

from snakemake.utils import R

rule eval_panel_global:
    input:  expand(DATA + 'interim/{dat}/{dat}.limit3.dat', dat=('clinvar', 'clinvar_single', 'clinvar_mult', 'clinvar_exp', 'denovo')),
            DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/uc.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/other/other.eff.dbnsfp.anno.hHack.dat.limit.xls'
    output: WORK + 'global.eval_panel.{cols}.stats',
            WORK + 'global.eval_panel.{cols}.eval',
            WORK + 'roc_df/{cols}'
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

def read_df(afile):
    col = afile.split('/')[-1]
    df_pre = pd.read_csv(afile, sep='\t')
    if col in df_pre.columns.values:
        df = df_pre[['chrom', 'pos', 'alt', col, 'y', 'Disease']]
    elif col + '_pred_lm' in df_pre.columns.values:
        df = df_pre[['chrom', 'pos', 'alt', col + '_pred_lm', 'y', 'Disease']]
    else:
        print(col)
    return df

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

feats = ('mpc', 'revel', 'ccr', 'is_domain')
combo_feats = ['-'.join(x) for x in powerset(feats) if x]
rule combine_predictions:
    input:  expand( WORK + 'roc_df/{cols}', cols=combo_feats )
    output: o = WORK + 'roc_df_combo'
    run:
        dfs = [read_df(afile) for afile in list(input)]
        keys = ['Disease', 'chrom', 'pos', 'alt', 'y']
        m = reduce(lambda left, right: pd.merge(left, right, on=keys), dfs)
        m.to_csv(output.o, index=False, sep='\t')

rule all_eval:
    input: expand( DOCS + 'plot/{method}.eval_panel.{cols}.totWrong.png', method=('global',), cols=combo_feats)

# ggplot(data=d) + geom_col(aes(y=var_count,x=score_type, fill=score_type)) + facet_grid(disease~., scale='free') + theme_bw() + ylab('Wrong Predictions') + theme(axis.text.x = element_text(angle=90)) + xlab('') + theme(legend.position="none")    


