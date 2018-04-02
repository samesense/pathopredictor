"""Pad benign variants.
   Train panel, and predict clinvar
"""

rule single_eval:
    input:  expand(DATA + 'interim/{dat}/{dat}.limit3.dat', dat=('clinvar', 'clinvar_single', 'clinvar_mult', 'clinvar_exp',)),
            DATA + 'interim/epi/EPIv6.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/epi/uc.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/other/other.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/gnomad/gnomad.rare.panel_{plow}_{phigh}',
    output:
            WORK + 'single_roc_df_clinvar/{cols}_{plow}_{phigh}'
    shell:  'python {SCRIPTS}single_gene.py {wildcards.cols} {input} {output}'

rule all_singles:
    input:  i = WORK + 'single_roc_df_clinvar/mpc-revel-ccr-is_domain_{plow}_{phigh}'
    output: o = WORK + 'single_{plow}_{phigh}.txt'
    run:
        diseases = {'genedx-epi-limitGene':'Epilepsy (dominant genes)',
                    'Rasopathies':'Rasopathies',
                    'genedx-epi':'Epilepsy',
                    'Cardiomyopathy':'Cardiomyopathy'}
        dd = pd.read_csv(input.i, sep='\t')
        crit = dd.apply(lambda row: 'tot' in row['Disease'] and not 'uc' in row['Disease'] and not 'ear' in row['Disease'], axis=1)
        df = dd[crit]
        df.loc[:, 'varType'] = df.apply(lambda row: 'Pathogenic' if 'ath' in row['PredictionStatusMPC'] else 'Benign', axis=1)
        df.loc[:, 'd'] = df.apply(lambda row: diseases[row['Disease'].split(':')[0]], axis=1)
        df.loc[:, 'key'] = df.apply(lambda row: row['d'] + ':' + row['gene'], axis=1)
        benign_genes = set(df[df.varType=='Benign']['key'])
        path_genes = set(df[df.varType=='Pathogenic']['key'])
        dat = []
        for key in path_genes - benign_genes:
            dis, gene = key.split(':')
            d = {'gene':gene, 'size':0, 'd':dis, 'PredictionStatusMPC':'CorrectBenign', 'varType':'Benign'}
            dat.append(d)
        ff = pd.DataFrame(dat)
        cols = ['d', 'PredictionStatusMPC', 'gene', 'varType']
        final = ( df[cols].groupby(cols)
                  .size().reset_index()
                  .rename(columns={0:'size'}) )
        df = pd.concat([final, ff])
        ( df.sort_values(by=['d', 'varType', 'gene', 'size'])
            .to_csv(output.o, index=False, sep='\t') )

rule plot_single_gene:
    input:  WORK + 'single_{plow}_{phigh}.txt'
    output: DOCS + 'paper_plots/single_gene_{plow}_{phigh}.pdf'
    run:
        R("""require(ggplot2)
             dat = read.delim("{input}", header=TRUE, sep="\t")
             p = ggplot(data=dat) + geom_col(aes(x=gene,y=size,fill=PredictionStatusMPC), position="dodge", stat="identity") +
             facet_wrap(d~varType, scale="free", ncol=2) + theme_bw(base_size=18) + coord_flip() +
             theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.y=element_blank(), axis.title.x=element_blank() )
             ggsave("{output}", p, width=20, height=20)""")

rule single:
    input: DOCS + 'paper_plots/single_gene_.003_.1.pdf', DOCS + 'paper_plots/single_gene_.003_1.pdf',


