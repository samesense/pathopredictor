"""Trained with panels,
   evaluate on clinvar
"""

rule eval_by_gene_clinvar:
    input:  i = WORK + 'roc_df_{eval_source}/{features}'
    output: o = DATA + 'interim/pred_clinvar_eval/{eval_source}.{features}.{evidenceCutoff}.{varTypes}'
    run:
        keys = ['Disease']
        df_pre = pd.read_csv(input.i, sep='\t')
        crit = df_pre.apply(lambda row: not 'uc' in row['Disease'] and not 'issue' in row['Disease'] and not 'earing' in row['Disease'], axis=1)
        df = df_pre[crit]

        apply_calc_wrong_baseline = mk_calc_wrong_func(int(wildcards.evidenceCutoff), 'PredictionStatusBaseline', wildcards.varTypes)
        apply_calc_wrong= mk_calc_wrong_func(int(wildcards.evidenceCutoff), 'PredictionStatusMPC', wildcards.varTypes)

        # baseline
        wrong_baseline_df = df.groupby(keys).apply(apply_calc_wrong_baseline).reset_index().rename(columns={0:'wrongFrac'})
        wrong_baseline_df['combo'] = 'clinvar_base.' + wildcards.features

        wrong_df= df.groupby(keys).apply(apply_calc_wrong).reset_index().rename(columns={0:'wrongFrac'})
        wrong_df['combo'] = 'clinvar.' + wildcards.features

        tot_wrong_df = df.groupby(['Disease']).apply(apply_calc_wrong).reset_index().rename(columns={0:'predictorWrongFracTot'})
        size_df = df.groupby(keys).size().reset_index().rename(columns={0:'size'})

        m0 = pd.concat([wrong_df, wrong_baseline_df])
        m = pd.merge(m0, size_df, on=keys)
        m.loc[:, 'clinvar_type'] = m.apply(lambda row: row['Disease'].split(':')[1], axis=1)
        m.loc[:, 'dd'] = m.apply(lambda row: row['Disease'].split(':')[0], axis=1)
        final_cols = ['Disease', 'dd', 'clinvar_type', 'combo',
                      'size', 'predictorWrongFracTot', 'wrongFrac']
        pd.merge(m, tot_wrong_df, on='Disease', how='left')[final_cols].to_csv(output.o, index=False, sep='\t')

def color_clinvar_bar(row):
    if row['combo'] in ('clinvar_base.mpc', 'clinvar_base.revel', 'clinvar_base.ccr', 'clinvar_base.is_domain'):
        return 'Baseline'
    if 'base' in row['combo']:
        return 'Combined baseline'
    return 'Trained'

def calc_clinvar_best_label(row):
    if row['panel_best'] and not row['is_best_clinvar']:
        return 'Best panel'
    if row['is_best_clinvar'] and not row['panel_best']:
        return 'Best ClinVar'
    if row['is_best_clinvar'] and row['panel_best']:
        return 'Best both'
    return 'NA'

rule combine_features_by_gene_clinvar_plot:
    input:  clinvar_data = expand(DATA + 'interim/pred_clinvar_eval/{{eval_source}}.{feature}.10.both', feature=COMBO_FEATS),
            panel_data = WORK + 'cc'
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
        dfp = pd.concat([pd.read_csv(afile, sep='\t') for afile in list(input.clinvar_data)])
        panel_cols = ['disease', 'st', 'is_best']
        panel_df = pd.read_csv(input.panel_data, sep='\t')[panel_cols].drop_duplicates().rename(columns={'disease':'dd', 'st':'combo','is_best':'panel_best'})
        crit = dfp.apply(lambda row: row['clinvar_type'] in clinvar_names, axis=1)
        df2 = dfp[crit]
        min_df = (df2[['clinvar_type', 'dd', 'wrongFrac']]
                  .groupby(['clinvar_type', 'dd'])
                  .apply(min).rename(columns={'dd':'dd_junk', 'clinvar_type':'clinvar_type_junk', 'wrongFrac':'min'})
                  .reset_index() )
        df3 = pd.merge(df2, min_df, on=['clinvar_type', 'dd'], how='left')
        df3.loc[:, 'Classifier'] = df3.apply(color_clinvar_bar, axis=1)
        df3.loc[:, 'combo'] = df3.apply(lambda row: row['combo'].replace('clinvar.', 'TRAINED_').replace('clinvar_base.', 'BASE_'), axis=1)
        df = pd.merge(df3, panel_df, on=['dd','combo'], how='left')
        df.loc[:, 'dd'] = df.apply(lambda row: diseases[row['dd']], axis=1)
        df.loc[:, 'clinvar_type'] = df.apply(lambda row: clinvar_names[row['clinvar_type']], axis=1)
        df.loc[:, 'is_best_clinvar'] = df.apply(lambda row: row['wrongFrac']==row['min'], axis=1)
        df.loc[:, 'is_best'] = df.apply(lambda row: row['panel_best'] or row['is_best_clinvar'], axis=1)
        df.loc[:, 'best_label'] = df.apply(calc_clinvar_best_label, axis=1)
        df.to_csv(output.o, index=False, sep='\t')

rule plot_clinvar_eval_paper:
    input:  i = DATA + 'interim/clinvar.by_gene_feat_combo.predictFullClinvar'
    output: DOCS + 'paper_plts/fig4_eval_clinvar.pdf'
    run:
        plot_cmd = """geom_col( aes(fill=Classifier, y=wrongFrac, x=reorder(combo, wrongFrac)) ) +
                      geom_point(data=dbest, aes(colour=best_label, x=combo,y=wrongFrac)) +
                      geom_text(data=label_df, aes(x=x,y=y,label=label))"""
        df = pd.read_csv(input.i, sep='\t')[['clinvar_type','dd','size']].drop_duplicates()
        df.loc[:, 'label'] = df.apply(lambda row: 'n=%d' % (row['size']), axis=1)
        df['y'] = 0.45
        df['x'] = 'TRAINED_mpc-revel-ccr'
        df.to_csv('tmp.clinvar.labels', index=False, sep='\t')

        R("""
          require(ggplot2)
          cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
          cbbPalette <- c("#D55E00", "#009E73", "#000000")
          d = read.delim("{input}", sep='\t', header=TRUE)
          label_df = read.delim("tmp.clinvar.labels", sep="\t", header=TRUE)
          dbest = d[d$is_best=="True",]
          d$clinvar_type = factor(d$clinvar_type, levels=c("Total ClinVar", "ClinVar w/ Evidence"))
          p = ggplot(data=d) + {plot_cmd} + scale_fill_manual(values=cbPalette) +
              ylab('Incorrect prediction fraction') + xlab('') + theme_bw(base_size=18) + facet_grid(clinvar_type~dd) +
              coord_flip() + labs(colour="") +
              theme(legend.position="bottom", axis.text.y = element_text(size=10)) + scale_colour_manual(values=cbbPalette)
          ggsave("{output}", p, height=9, width=20)
          """)
        shell('rm tmp.clinvar.labels')

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
    input: 
           DOCS + 'plot/eval_clinvar/clinvar.genedx-epi:tot.predictFullClinvar.byVarClassFalse.png'

rule all_clinvar_eval:
    input:
           expand(DOCS + 'plot/eval_clinvar/{eval_source}.{disease}:{d2}.predictFullClinvar.png', eval_source=('clinvar',), d2=('tot', ), disease=DD), \



