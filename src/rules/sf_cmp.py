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
            q = DATA + 'interim/enrich_q/{var}.xls'
    output: o = DATA + 'interim/sig_flagged/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.{var}.dat'
    run:
        use_cols = ['pfam', 'fg_gtr', 'qval']
        r = {'fg_gtr':wildcards.var + '_fg_gtr',
             'qval':wildcards.var + '_qval'}
        q_df = pandas.read_csv(input.q, sep='\t', usecols=use_cols).rename(columns=r)
        df = pandas.read_csv(input.d, sep='\t')
        m = pandas.merge(df, q_df, on='pfam', how='left').fillna(0)
        m.to_csv(output.o, index=False, sep='\t')

rule test_pathogenic_enrichment:
    input:  DATA + 'interim/sig_flagged/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.{var}.dat'
    output: DATA + 'interim/path_enrich/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.{var}.dat'
    shell:  'python {SCRIPTS}eval_pathogenic_enrichment.py {input} {wildcards.var} {output}'

vars = ('disruptive_inframe_insertion',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'disruptive_inframe_deletion',
        'frameshift_variant',
        'stop_gained',
        'frameshift_variant+start_lost',
        'initiator_codon_variant',
        'stop_lost+inframe_deletion',
        'frameshift_variant+stop_lost',
        'stop_gained+disruptive_inframe_deletion',
        'stop_lost',
        'frameshift_variant+stop_gained',
        'stop_retained_variant',
        'start_lost',
        'splice_region_variant',
        'mis', 'all', 'lof')

rule combine_path_burden:
    input: expand(DATA + 'interim/path_enrich/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.{var}.dat', \
                  var=vars)
    output: DATA + 'interim/path_enrich_eval'
    run:
        f = list(input)[0]
        shell('head -1 {f} > {output}')
        for l in list(input):
            shell('grep -v benign {l} >> {output}')

rule plot:
    input:  DATA + 'interim/path_enrich_eval'
    output: PLOTS + 'path_frac_wo_vus.png',
            PLOTS + 'path_frac_w_vus.png',
            PLOTS + 'benign_frac_w_vus.png',
            PLOTS + 'var_counts_by_enrichment.png'
    run:  
        shell('Rscript {SCRIPTS}plot_fracs.R {input} {output}')
        shell('rm Rplots.pdf')

rule all_cmp:
    input: DATA + 'interim/sig_flagged/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.all.dat', PLOTS + 'var_counts_by_enrichment.png'
