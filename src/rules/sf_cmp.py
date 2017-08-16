"""Fisher tests for burdened domains"""
import pandas
import statsmodels.sandbox.stats.multicomp as fdr
from const import *

include: "sf_exac.py"
include: "sf_ann.py"

rule merge:
    input:  DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.dat',
            DATA + 'interim/r1_no_tcga/exac.tidy.eff.dbnsfp.gt.anno.hHack.splitPfam.dat'
    output: DATA + 'interim/merge/{var}'
    shell:  'python {SCRIPTS}merge.py {wildcards.var} {input} {output}'

rule test:
    input:  DATA + 'interim/merge/{var}'
    output: DATA + 'interim/enrich/{var}'
    shell:  '{PY27_T} {SCRIPTS}perry_fisher.py {input} {output}'

rule fdr:
    input:  d = DATA + 'interim/enrich/{var}'
    output: o = DATA + 'interim/enrich_q/{var}.xls'
    run:
        df = pandas.read_csv(input.d, sep='\t')
        (reject_ls, qs, alpha_sidak, alpha_bond) = fdr.multipletests(df['pval'].values, method='fdr_bh')
        df['qval'] = qs
        df.sort_values(by='qval', ascending=True).to_csv(output.o, index=False, sep='\t')

rule flag_sig_domains:
    input:  d = DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.dat',
            q = DATA + 'interim/enrich_q/missense_variant.xls'
    output: o = DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.missense.dat'
    run:
        use_cols = ['pfam', 'fg_gtr', 'qval']
        r = {'fg_gtr':'missense_fg_gtr',
             'qval':'missense_qval'}
        q_df = pandas.read_csv(input.q, sep='\t', usecols=use_cols).rename(columns=r)
        df = pandas.read_csv(input.d, sep='\t')
        m = pandas.merge(df, q_df, on='pfam', how='left')
        m.to_csv(output.o, index=False, sep='\t')

rule m:
    input: expand(DATA + 'interim/enrich_q/{var}.xls', var=('missense_variant', 'stop_gained'))
