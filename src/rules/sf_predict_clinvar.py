"""Trained with panels,
   evaluate on clinvar
"""

rule eval_by_gene_clinvar:
    input:  i = WORK + 'roc_df_{eval_source}/{features}'
    output: o = DATA + 'interim/pred_clinvar_eval/{eval_source}.{features}'
    run:
        keys = ['Disease']
        df_pre = pd.read_csv(input.i, sep='\t')
        crit = df_pre.apply(lambda row: not 'uc' in row['Disease'] and not 'issue' in row['Disease'] and not 'earing' in row['Disease'], axis=1)
        df = df_pre[crit]
        wrong_df= df.groupby(keys).apply(calc_wrong).reset_index().rename(columns={0:'wrongFrac'})
        tot_wrong_df = df.groupby(['Disease']).apply(calc_wrong).reset_index().rename(columns={0:'predictorWrongFracTot'})
        size_df = df.groupby(keys).size().reset_index().rename(columns={0:'size'})
        m = pd.merge(wrong_df, size_df, on=keys)
        m.loc[:, 'clinvar_type'] = m.apply(lambda row: row['Disease'].split(':')[1], axis=1)
        m.loc[:, 'dd'] = m.apply(lambda row: row['Disease'].split(':')[0], axis=1)
        pd.merge(m, tot_wrong_df, on='Disease', how='left')[['Disease', 'dd', 'clinvar_type', 'size', 'predictorWrongFracTot', 'wrongFrac']].to_csv(output.o, index=False, sep='\t')

def color_clinvar_bar(row):
    if row['combo'] in ('paper_mpc', 'paper_revel', 'paper_ccr'):
        return 'Baseline'
    if 'paper' in row['combo']:
        return 'Combined baseline'
    return 'Trained'

rule combine_features_by_gene_clinvar_plot:
    input:  expand(DATA + 'interim/pred_clinvar_eval/{{eval_source}}.{feature}', feature=COMBO_FEATS)
    output: o=DATA + 'interim/{eval_source}.by_gene_feat_combo.predictFullClinvar'
    run:
        disease_order = {'genedx-epi-limitGene':3,
                         'Rasopathies':4,
                         'genedx-epi':2,
                         'Cardiomyopathy':1}

        diseases = {'genedx-epi-limitGene':'Epilepsy (dominant genes)',
                    'Rasopathies':'Rasopathies',
                    'genedx-epi':'Epilepsy',
                    'Cardiomyopathy':'Cardiomyopathy'}
        clinvar_names  = {'tot': 'Total ClinVar',
                          'single': 'ClinVar w/ Evidence'}
        dfp = pd.concat([read_gene_df(afile) for afile in list(input)])
        crit = dfp.apply(lambda row: row['clinvar_type'] in clinvar_names, axis=1)
        df2 = dfp[crit]
        min_df = (df2[['clinvar_type', 'dd', 'wrongFrac']]
                  .groupby(['clinvar_type', 'dd'])
                  .apply(min).rename(columns={'dd':'dd_junk', 'clinvar_type':'clinvar_type_junk', 'wrongFrac':'min'})
                  .reset_index() )
        df = pd.merge(df2, min_df, on=['clinvar_type', 'dd'], how='left')
        df.loc[:, 'dd'] = df.apply(lambda row: diseases[row['dd']], axis=1)
        df.loc[:, 'clinvar_type'] = df.apply(lambda row: clinvar_names[row['clinvar_type']], axis=1)
        df.loc[:, 'Classifier'] = df.apply(color_clinvar_bar, axis=1)
        df.loc[:, 'combo'] = df.apply(lambda row: row['combo'].replace('clinvar.', 'TRAINED_'), axis=1)
        df.loc[:, 'is_best'] = df.apply(lambda row: row['wrongFrac']==row['min'], axis=1)
        df.to_csv(output.o, index=False, sep='\t')

rule plot_clinvar_eval_paper:
    input:  DATA + 'interim/clinvar.by_gene_feat_combo.predictFullClinvar'
    output: DOCS + 'paper_plts/fig4_eval_clinvar.pdf'
    run:
        plot_cmd = 'geom_col( fill="#56B4E9", aes(y=wrongFrac, x=reorder(combo, predictorWrongFracTot)) ) + geom_point(data=dbest, aes(x=combo,y=wrongFrac))'
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          dbest = d[d$is_best=="True",]
          d$clinvar_type = factor(d$clinvar_type, levels=c("Total ClinVar", "ClinVar w/ Evidence"))
          p = ggplot(data=d) + {plot_cmd} +
              ylab('Incorrect prediction fraction') + xlab('') + theme_bw() + facet_grid(clinvar_type~dd) +
              coord_flip() +
              theme(axis.text.y = element_text(size=10))
          ggsave("{output}", p, height=10, width=20)
          """)

rule plot_clinvar_eval:
    input:  DATA + 'interim/{eval_source}.by_gene_feat_combo.predictFullClinvar'
    output: DOCS + 'plot/eval_clinvar/{eval_source}.{disease}.predictFullClinvar.byVarClass{byVarClass}.png'
    run:
        dd = wildcards.disease.split(':')[0]
        if wildcards.byVarClass == 'True':
            plot_cmd = 'geom_col(aes(y=wrongFrac, x=reorder(combo, predictorWrongFracTot), fill=var_class), position="dodge")'
        else:
            plot_cmd = 'geom_col(aes(y=wrongFrac, x=reorder(combo, predictorWrongFracTot)), position="dodge")'
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          p = ggplot(data=d[d$dd=="{dd}",]) + {plot_cmd} +
              ylab('Wrong prediction fraction') + xlab('') + theme_bw() + facet_grid(clinvar_type~.) +
              labs(fill="Variant class") + coord_flip() +
              theme(axis.text.y = element_text(size=6))
          ggsave("{output}", p, height=20)
          """)

rule tmp_c:
    input: DOCS + 'plot/eval_clinvar/clinvar.genedx-epi:tot.predictFullClinvar.byVarClassTrue.png',
           DOCS + 'plot/eval_clinvar/clinvar.genedx-epi:tot.predictFullClinvar.byVarClassFalse.png'

rule all_clinvar_eval:
    input:
           expand(DOCS + 'plot/eval_clinvar/{eval_source}.{disease}:{d2}.predictFullClinvar.png', eval_source=('clinvar',), d2=('tot', ), disease=DD), \



