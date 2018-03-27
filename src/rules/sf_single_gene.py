"""Pad benign variants.
   Train panel, and predict clinvar
"""

rule single_eval:
    input:  expand(DATA + 'interim/{dat}/{dat}.limit3.dat', dat=('clinvar', 'clinvar_single', 'clinvar_mult', 'clinvar_exp',)),
            DATA + 'interim/epi/EPIv6.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/epi/uc.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/other/other.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/gnomad/gnomad.rare.panel',
            DATA + 'interim/gnomad/gnomad.rare.clinvar'
    output:
            WORK + 'single_roc_df_clinvar/{cols}'
    shell:  'python {SCRIPTS}single_gene.py {wildcards.cols} {input} {output}'

rule all_singles:
    input:  i = WORK + 'single_roc_df_clinvar/mpc-revel-ccr-is_domain'
    output: o = WORK + 'single.txt'
    run:
        df = pd.read_csv(input.i, sep='\t')
        crit = df.apply(lambda row: 'tot' in row['Disease'] and not 'uc' in row['Disease'] and not 'ear' in row['Disease'], axis=1)
        cols = ['Disease', 'PredictionStatusMPC', 'gene']
        df[crit][cols].groupby(cols).size().reset_index().rename(columns={0:'size'}).sort_values(by=['Disease', 'gene', 'size']).to_csv(output.o, index=False, sep='\t')
