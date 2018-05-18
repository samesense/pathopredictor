"""Predictor paper plots."""

rule count_plot_data:
    input: panel = WORK + 'roc_df_panel/' + C_FEATS,
           clinvar = WORK + 'roc_df_clinvar/' + C_FEATS
    output: o = WORK + 'paper_plot_data/count_plot'
    run:
        diseases = {'genedx-epi-limitGene':'Epilepsy (dominant genes)',
                    'Rasopathies':'Rasopathies',
                    'EPI':'Epilepsy',
                    'Cardiomyopathy':'Cardiomyopathy'}
        clin_labels = {'single':'ClinVar w/ Evidence',
                       'tot': 'Total ClinVar'}
        df_panel = pd.read_csv(input.panel, sep='\t')
        panel = df_panel.groupby(['Disease','y']).size().reset_index().rename(columns={0:'count'})
        panel['eval_type'] = 'Panel'
        panel.loc[:, 'Disease'] = panel.apply(lambda row: diseases[row['Disease']], axis=1)

        df_clinvar = pd.read_csv(input.clinvar, sep='\t')
        df_clinvar_single = df_clinvar[df_clinvar.is_single]

        clinvar = df_clinvar.groupby(['Disease','y']).size().reset_index().rename(columns={0:'count'})
        clinvar.loc[:, 'eval_type'] = clin_labels['tot']
        clinvar.loc[:, 'Disease'] = clinvar.apply(lambda row: diseases[row['Disease']], axis=1)

        s_clinvar = df_clinvar_single.groupby(['Disease','y']).size().reset_index().rename(columns={0:'count'})
        s_clinvar.loc[:, 'eval_type'] = clin_labels['single']
        s_clinvar.loc[:, 'Disease'] = s_clinvar.apply(lambda row: diseases[row['Disease']], axis=1)
        
        df = pd.concat([panel, clinvar, s_clinvar])
        df.loc[:,'y'] = df.apply(lambda row: 'Pathogenic' if row['y']==1 else 'Benign', axis=1)
        df['count_type'] = 'Variants'

        panel_gene = df_panel[['Disease', 'gene']].drop_duplicates().groupby(['Disease']).size().reset_index().rename(columns={0:'count'})
        panel_gene['eval_type'] = 'Panel'
        panel_gene.loc[:, 'Disease'] = panel_gene.apply(lambda row: diseases[row['Disease']], axis=1)
        
        clinvar_gene = df_clinvar[['Disease','gene']].drop_duplicates().groupby(['Disease']).size().reset_index().rename(columns={0:'count'})
        clinvar_gene.loc[:, 'eval_type'] = clin_labels['tot']
        clinvar_gene.loc[:, 'Disease'] = clinvar_gene.apply(lambda row: diseases[row['Disease']], axis=1)

        s_clinvar_gene = df_clinvar_single[['Disease','gene']].drop_duplicates().groupby(['Disease']).size().reset_index().rename(columns={0:'count'})
        s_clinvar_gene.loc[:, 'eval_type'] = clin_labels['single']
        s_clinvar_gene.loc[:, 'Disease'] = s_clinvar_gene.apply(lambda row: diseases[row['Disease']], axis=1)
        
        df_gene = pd.concat([panel_gene, clinvar_gene, s_clinvar_gene])
        df_gene['y'] = 'Gene'
        df_gene['count_type'] = 'Genes'
        pd.concat([df, df_gene]).to_csv(output.o, index=False, sep='\t')

rule count_plot:
    input:  WORK + 'paper_plot_data/count_plot'
    output: DOCS + 'paper_plts/fig1_countPlot.pdf'
    run:
        R("""
          require(ggplot2)
          require(forcats)
          d = read.delim("{input}", sep='\t', header=TRUE)
          d$eval_type = factor(d$eval_type, levels=c("Panel", "Total ClinVar", "ClinVar w/ Evidence"))
          p = ggplot(data=d) +
              geom_bar(stat="identity", aes(x=Disease, y=count, fill=fct_reorder(y,count))) +
              facet_grid(count_type~eval_type, scale="free_y") + theme_bw(base_size=18) +
              ylab('Count') + labs(fill= "") + theme(legend.position="bottom") +
              theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank())
          ggsave("{output}", p, width=10)
          """)

rule combine_heatmap_clinvar_and_panel:
    input:  clinvar = DATA + 'interim/clinvar.by_gene_feat_combo.{evidenceCutoff}.{varTypes}',
            panel = DATA + 'interim/panel.by_gene_feat_combo.{evidenceCutoff}.{varTypes}'
    output: o = WORK + 'paper_plot_data/heatmap.{evidenceCutoff}.{varTypes}'
    run:
        diseases = {'genedx-epi-limitGene':'Epilepsy (dominant genes)',
                    'Rasopathies':'Rasopathies',
                    'genedx-epi':'Epilepsy',
                    'Cardiomyopathy':'Cardiomyopathy'}
        panel_df_pre = pd.read_csv(input.panel, sep='\t')
        clinvar_df_pre = pd.read_csv(input.clinvar, sep='\t')
        crit = panel_df_pre.apply(lambda row: row['Disease'] in diseases, axis=1)
        panel_df = panel_df_pre[crit]
        panel_df.loc[:, 'disease'] = panel_df.apply(lambda row: diseases[row['Disease']],
                                                    axis=1)
        crit = clinvar_df_pre.apply(lambda row: row['Disease'].split(':')[0]
                                    in diseases, axis=1)
        clinvar_df = clinvar_df_pre[crit]
        clinvar_df.loc[:, 'disease'] = clinvar_df.apply(lambda row:
                                                        diseases[row['Disease'].split(':')[0]], axis=1)
        df = ( pd.concat([clinvar_df[['disease','gene']].drop_duplicates(),
                          panel_df[['disease','gene']].drop_duplicates()])
                 .groupby(['disease', 'gene']).size().reset_index()
                 .rename(columns={0:'size'})
             )
        disease_gene_keep = {row['disease'] + ':' + row['gene']
                             for _, row in df[df['size']==2].iterrows()}
        crit = panel_df.apply(lambda row: row['disease'] + ':' + row['gene']
                              in disease_gene_keep, axis=1)
        panel_final = panel_df[crit][['disease', 'combo', 'gene', 'wrongFrac', 'var_path_count', 'var_benign_count']]
        panel_final['panel_or_clinvar'] = 'Panel'
        panel_final.loc[:, 'gene'] = panel_final.apply(lambda row: row['gene'] + ' Panel', axis=1)
        min_df = (panel_final[['disease', 'gene', 'wrongFrac']]
                  .groupby(['disease', 'gene']).apply(min)
                  .rename(columns={'gene':'gene_junk', 'disease':'dis_junk', 'wrongFrac':'min'})
                  .reset_index())
        pm = pd.merge(panel_final, min_df, on=['disease', 'gene'])
        pm.loc[:, 'is_best'] = pm.apply(lambda row: row['wrongFrac'] == row['min'], axis=1)

        crit = clinvar_df.apply(lambda row: row['disease'] + ':' + row['gene']
                                in disease_gene_keep, axis=1)
        clinvar_final = clinvar_df[crit][['disease', 'combo', 'gene', 'wrongFrac', 'var_path_count', 'var_benign_count']]
        clinvar_final['panel_or_clinvar'] = 'Total ClinVar'
        clinvar_final.loc[:, 'gene'] = clinvar_final.apply(lambda row: row['gene'] + ' Total ClinVar', axis=1)
        min_df = ( clinvar_final[['disease', 'gene', 'wrongFrac']]
                   .groupby(['disease', 'gene']).apply(min)
                   .rename(columns={'gene':'gene_junk', 'disease':'dis_junk', 'wrongFrac':'min'})
                   .reset_index() )
        pc = pd.merge(clinvar_final, min_df, on=['disease', 'gene'])
        pc.loc[:, 'is_best'] = pc.apply(lambda row: row['wrongFrac'] == row['min'], axis=1)

        df_final = pd.concat([pc, pm])
        df_final.loc[:, 'combo'] = df_final.apply(lambda row: row['combo'].replace('_base.', 'BASE ').replace('_trained.','TRAINED '), axis=1)
        df_final.to_csv(output.o, index=False, sep='\t')

rule paper_heatmap:
    input:  WORK + 'paper_plot_data/heatmap.{evidenceCutoff}.{varTypes}'
    output: DOCS + 'paper_plts/fig5_heatmap.{evidenceCutoff}.{varTypes}.pdf'
    run:
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          best_d = d[d$is_best=="True",]
          p = ggplot(data=d) +
              geom_raster(aes(y=gene, x=reorder(combo, wrongFrac, mean), fill=wrongFrac)) +
              geom_point(data=best_d, aes(y=gene, x=combo)) +
              ylab('') + xlab('') + theme_bw(base_size=18) +
              scale_fill_gradient(low="yellow", high="blue") +
              labs(fill="Incorrect prediction fraction") +
              facet_grid(disease~., scale="free") +
              theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
          ggsave("{output}", p, width=20, height=20)
          """)

rule mk_size_plot_data:
    input:  i = WORK + 'paper_plot_data/heatmap.{evidenceCutoff}.{varTypes}'
    output: o = WORK + 'paper_plot_data/size.{evidenceCutoff}.{varTypes}'
    run:
        df = pd.read_csv(input.i, sep='\t')
        cols = set(df.columns.values)
        path_cols = list(cols - set(('var_benign_count')))
        benign_cols = list(cols - set(('var_path_count')))
        path_df = df[path_cols].rename(columns={'var_path_count':'var_count'})
        path_df['VariantType'] = 'Pathogenic'

        benign_df = df[benign_cols].rename(columns={'var_benign_count':'var_count'})
        benign_df['VariantType'] = 'Benign'

        df = pd.concat([path_df, benign_df])
        df[['VariantType', 'gene', 'disease', 'var_count']].drop_duplicates().to_csv(output.o, index=False, sep='\t')

rule size_bar_paper_plot:
    input:  WORK + 'paper_plot_data/size.{evidenceCutoff}.{varTypes}'
    output: DOCS + 'paper_plts/fig5b.varCount.{evidenceCutoff}.{varTypes}.pdf'
    run:
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          p = ggplot(data=d) +
              geom_col(aes(fill=VariantType, x=gene, y=var_count), stat="identity") + facet_grid(disease~., scale="free") +
              coord_flip() + xlab('') + ylab('Variant count') + theme_bw(base_size=18) + labs(fill='')
          ggsave("{output}", p, height=20)
          """)

#FIGS = ('fig1_count_plot', 'fig4_eval_clinvar', 'fig3_panelEval.byVarClassFalse', 'fig5_idi', 'fig6_single_gene_collapse_.003_.1')
#FIGS = ('fig1_count_plot', 'fig3_panelEval.byVarClassFalse')
FIGS = ('fig1_countPlot', 'fig5_evalClinvar', 'fig4_panelEval',
        'fig3b_featureCor', 'fig3a_featureImportance',
        'fig6_single_gene_collapse_.003_1.full',
        'fig7_single_gene_collapse_.003_1.single', )

rule upload_paper_plot:
    input:  DOCS + 'paper_plts/{fig}.pdf'
    output: DBox.remote('ahmad_predictor/{fig}.pdf')
    shell:  'cp {input} {output}'

rule all_paper_plots:
    input: expand(DOCS + 'paper_plts/{fig}.pdf', fig=FIGS)

rule upload_all:
    input: expand(DBox.remote('ahmad_predictor/{fig}.pdf'), fig=FIGS)
           #expand(DOCS + 'paper_plts/fig5_heatmap.{evidenceCutoff}.{varTypes}.pdf', varTypes=('pathogenic', 'benign', 'both', 'bothAhmad'), evidenceCutoff=(4,5,10)),
           #expand(DOCS + 'paper_plts/fig5b.varCount.{evidenceCutoff}.{varTypes}.pdf', varTypes=('pathogenic', 'benign', 'both', 'bothAhmad'), evidenceCutoff=(4,5,10))
