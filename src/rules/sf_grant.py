"""Grant figs"""
import pandas

rule fg_lolly:
    input:  i = DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.xls'
    output: DOCS + 'plots/{gene}.panel_one.lolly.png'
    run:  
        benign_color = '00bfc4'
        path_color = 'f8766d'
        df_pre = pandas.read_csv(input.i, sep='\t')
        df = df_pre[ (df_pre.gene==wildcards.gene) & (df_pre.eff=='missense_variant') ]
        benign_ls = [x + '#' + benign_color for x in 
                     set(df[(df.clin_class=='BENIGN') | (df.clin_class=='LIKELY_BENIGN')]['Protein_Change'].values)]
        path_ls = [x + '#' + path_color for x in 
                   set(df[(df.clin_class=='PATHOGENIC') | (df.clin_class=='LIKLEY_PATHOGENIC')]['Protein_Change'].values)]
        shell('~/me/bin/lollipops -o={output} -f=/home/evansj/me/fonts/arial.ttf {wildcards.gene} {path_ls} {benign_ls}')

rule lollies:
    input: expand( DOCS + 'plots/{gene}.panel_one.lolly.png', gene=('SCN1A', 'TSC2', 'SPTAN1', 'TSC1', 'GRIN2A', 'PNPO', 'PPT1') )
