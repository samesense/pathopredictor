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
        (reject_ls, qs, alpha_sidak, alpha_bond) = fdr.multipletests(df['rare_pval'].values, method='fdr_bh')
        df['rare_qval'] = qs
        (reject_ls, qs, alpha_sidak, alpha_bond) = fdr.multipletests(df['rarem_pval'].values, method='fdr_bh')
        df['rarem_qval'] = qs
        df.sort_values(by='rare_qval', ascending=True).to_csv(output.o, index=False, sep='\t')

rule flag_sig_domains:
    input:  d = DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.dat',
            q = DATA + 'interim/enrich_q/{var}.xls'
    output: o = DATA + 'interim/sig_flagged/{var}/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.dat'
    run:
        v = ('fg_gtr', 'qval', 'pos_fam', 'fg_tot', 'ac', 'bg_tot')
        use_cols = ['pfam'] + ['rare_' + x for x in v] + ['rarem_' + x for x in v]
        r = {}
        for l in ('rare_', 'rarem_'):
            for vv in v:
                r[l + vv] = l + wildcards.var + '_' + vv
        # r = {'fg_gtr':wildcards.var + '_fg_gtr',
        #      'qval':wildcards.var + '_qval',
        #      'pos_fam':wildcards.var + '_fg',
        #      'fg_tot':wildcards.var + '_fgtot',
        #      'ac':wildcards.var + '_bg', 
        #      'bg_tot':wildcards.var + '_bgtot'}
        q_df = pandas.read_csv(input.q, sep='\t', usecols=use_cols).rename(columns=r)
        df = pandas.read_csv(input.d, sep='\t')
        m = pandas.merge(df, q_df, on='pfam', how='left').fillna(0)
        m.to_csv(output.o, index=False, sep='\t')

rule test_pathogenic_enrichment:
    input:  DATA + 'interim/sig_flagged/{eff}/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.dat'
    output: DATA + 'interim/path_enrich/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.limit_{limit}.eff_{eff}.dat'
    shell:  'python {SCRIPTS}eval_pathogenic_enrichment.py {input} {wildcards.eff} {wildcards.limit} {output}'

rule assign_sig:
    input:  DATA + 'interim/sig_flagged/{eff}/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.dat'
    output: DATA + 'interim/sig_flagged_assigned/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.limit_{limit}.eff_{eff}.dat'
    shell:  'python {SCRIPTS}eval_mpc_enrichment.py {input} {wildcards.eff} {wildcards.limit} {output}'

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

vars = ('mis', 'all',)

rule combine_path_burden:
    input: expand(DATA + 'interim/path_enrich/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.limit_{{limit}}.eff_{eff}.dat', \
                  eff=vars)
    output: DATA + 'interim/path_enrich_eval.{limit}'
    run:
        f = list(input)[0]
        shell('head -1 {f} > {output}')
        for l in list(input):
            shell('grep -v benign {l} >> {output}')

rule plot:
    input:  DATA + 'interim/path_enrich_eval.{limit}'
    output: PLOTS + '{limit}.path_frac_wo_vus.png',
            PLOTS + '{limit}.path_frac_w_vus.png',
            PLOTS + '{limit}.benign_frac_w_vus.png',
            PLOTS + '{limit}.var_counts_by_enrichment.png'
    run:  
        shell('Rscript {SCRIPTS}plot_fracs.R {input} {output}')
        shell('rm Rplots.pdf')

# DATA + 'interim/sig_flagged/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.all.dat',
rule all_cmp:
    input: PLOTS + 'rare.var_counts_by_enrichment.png',
           PLOTS + 'rarem.var_counts_by_enrichment.png',
           expand(DATA + 'interim/sig_flagged_assigned/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.limit_{limit}.eff_{var}.dat', \
                  limit=('rarem', 'rare'), var=vars)
