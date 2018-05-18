"""Pad benign variants.
   Train panel, and predict clinvar
"""

rule mk_single_gene_data:
    input:  DATA + 'interim/panel.dat',
            DATA + 'interim/clinvar.dat',
            DATA + 'interim/gnomad/gnomad.rare.panel_{plow}_{phigh}',
    output: DATA + 'interim/single_gene/dat_{plow}_{phigh}'
    shell:  'python {SCRIPTS}mk_single_gene_data.py {input} {output}'

rule tmp_single:
    input: DATA + 'interim/single_gene/dat_.003_1'

rule single_eval:
    input:  DATA + 'interm/single_gene/dat_{plow}_{phigh}'
    output: DATA + 'interim/single_gene/roc_df_clinvar/{cols}_{plow}_{phigh}'
    shell:  'python {SCRIPTS}single_gene.py {wildcards.cols} {input} {output}'

def score_gene(rows):
    disease = list(set(rows['d']))[0]
    pathogenic = len(rows[rows.varType=='Pathogenic'])
    benign = len(rows[rows.varType=='Benign'])
    label = 'p=%d,b=%d' % (pathogenic, benign)
    correct = len(rows[ (rows.PredictionStatusMPC=='CorrectBenign') | (rows.PredictionStatusMPC=='CorrectPath') ])
    accuracy = float(correct)/(benign + pathogenic)
    dat = {'disease':disease, 'count_label':label, 'accuracy':accuracy}
    s = pd.Series(dat, index=['disease', 'accuracy', 'count_label', 'accuracy'])
    return s

rule score_single_gene:
    input:  i = WORK + 'single_gene/roc_df_clinvar/{plow}_{phigh}/{features}'
    output: o = DATA + 'interim/single_gene/evals/{features}.single_{plow}_{phigh}.txt'
    run:
        diseases = {'genedx-epi-limitGene':'Epilepsy (dominant genes)',
                    'Rasopathies':'Rasopathies',
                    'genedx-epi':'Epilepsy',
                    'Cardiomyopathy':'Cardiomyopathy'}
        dd = pd.read_csv(input.i, sep='\t')
        df.loc[:, 'varType'] = df.apply(lambda row: 'Pathogenic' if 'ath' in row['PredictionStatusMPC'] else 'Benign', axis=1)
        df.loc[:, 'd'] = df.apply(lambda row: diseases[row['Disease'].split(':')[0]], axis=1)
        eval_df = df.groupby(['clinvar_set', 'gene']).apply(score_gene).reset_index()
        eval_df[:, 'terms'] = wildcards.features
        ( eval_df.sort_values(by=['accuracy',])
          .to_csv(output.o, index=False, sep='\t') )

rule collapse_single_gene:
    input:  expand(DATA + 'interim/single_gene/evals/{features}.single_{{plow}}_{{phigh}}.txt', features = COMBO_FEATS )
    output: DATA + 'interim/single_gene/plot_data_{plow}_{phigh}'
    run:
        pd.concat([pd.read_csv(afile, sep='\t') for afile in input]).to_csv(output.o, index=False, sep='\t')

rule plot_single_gene:
    input:  DATA + 'interim/single_gene/plot_data_{plow}_{phigh}'
    output: DOCS + 'paper_plts/fig6_single_gene_collapse_{plow}_{phigh}.pdf'
    run:
        R("""require(ggplot2)
             dat = read.delim("{input}", header=TRUE, sep="\t")
             p = ggplot(data=dat) + geom_col(aes(x=gene,y=accuracy,fill=terms), stat="identity", position="dodge") +
             facet_grid(d~clinvar_set, scale="free") + theme_bw(base_size=18) + coord_flip() +
             theme(axis.title.y=element_blank(), axis.title.x=element_blank() )
             ggsave("{output}", p, width=10, height=15)""")

rule single:
    input: DOCS + 'paper_plts/fig6_single_gene_collapse_.003_.1.pdf', DOCS + 'paper_plts/fig6_single_gene_collapse_.003_1.pdf',


