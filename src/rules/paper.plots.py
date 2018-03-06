"""Predictor paper plots."""

rule count_plot_data:
    input: panel = WORK + 'roc_df_panel/mpc',
           clinvar = WORK + 'roc_df_clinvar/mpc'
    output: o = WORK + 'paper_plot_data/count_plot'
    run:
        diseases = {'genedx-epi-limitGene':'Epilepsy (dominant genes)',
                    'Rasopathies':'Rasopathies',
                    'genedx-epi':'Epilepsy',
                    'Cardiomyopathy':'Cardiomyopathy'}
        clin_labels = {'single':'ClinVar w/ Evidence',
                       'tot': 'Total ClinVar'}
        df_panel = pd.read_csv(input.panel, sep='\t')
        crit_p = df_panel.apply(lambda row: row['Disease'] in diseases, axis=1)
        panel = df_panel[crit_p].groupby(['Disease','y']).size().reset_index().rename(columns={0:'count'})
        panel['eval_type'] = 'Panel'
        panel.loc[:, 'Disease'] = panel.apply(lambda row: diseases[row['Disease']], axis=1)
        df_clinvar = pd.read_csv(input.clinvar, sep='\t')
        crit_c = df_clinvar.apply(lambda row: row['Disease'].split(':')[0] in diseases and row['Disease'].split(':')[1] in ('tot', 'single'), axis=1)
        clinvar = df_clinvar[crit_c].groupby(['Disease','y']).size().reset_index().rename(columns={0:'count'})
        clinvar.loc[:, 'eval_type'] = clinvar.apply(lambda row: clin_labels[row['Disease'].split(':')[1]], axis=1)
        clinvar.loc[:, 'Disease'] = clinvar.apply(lambda row: diseases[row['Disease'].split(':')[0]], axis=1)
        df = pd.concat([panel, clinvar])
        df.loc[:,'y'] = df.apply(lambda row: 'Pathogenic' if row['y']==1 else 'Benign', axis=1)
        df['count_type'] = 'Variants'

        panel_gene = df_panel[crit_p][['Disease', 'gene']].drop_duplicates().groupby(['Disease']).size().reset_index().rename(columns={0:'count'})
        panel_gene['eval_type'] = 'Panel'
        panel_gene.loc[:, 'Disease'] = panel_gene.apply(lambda row: diseases[row['Disease']], axis=1)
        clinvar_gene = df_clinvar[crit_c][['Disease','gene']].drop_duplicates().groupby(['Disease']).size().reset_index().rename(columns={0:'count'})
        clinvar_gene.loc[:, 'eval_type'] = clinvar_gene.apply(lambda row: clin_labels[row['Disease'].split(':')[1]], axis=1)
        clinvar_gene.loc[:, 'Disease'] = clinvar_gene.apply(lambda row: diseases[row['Disease'].split(':')[0]], axis=1)
        df_gene = pd.concat([panel_gene, clinvar_gene])
        df_gene['y'] = 'Gene'
        df_gene['count_type'] = 'Genes'
        pd.concat([df, df_gene]).to_csv(output.o, index=False, sep='\t')

rule count_plot:
    input:  WORK + 'paper_plot_data/count_plot'
    output: DOCS + 'paper_plts/fig1_count_plot.pdf'
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
    input:  clinvar = DATA + 'interim/clinvar.by_gene_feat_combo.5',
            panel = DATA + 'interim/panel.by_gene_feat_combo.5'
    output: o = WORK + 'paper_plot_data/heatmap'
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
        panel_final = panel_df[crit][['disease', 'combo', 'gene', 'wrongFrac']]
        panel_final['panel_or_clinvar'] = 'Panel'
        panel_final.loc[:, 'gene'] = panel_final.apply(lambda row: row['gene'] + ' Panel', axis=1)

        crit = clinvar_df.apply(lambda row: row['disease'] + ':' + row['gene']
                                in disease_gene_keep, axis=1)
        clinvar_final = clinvar_df[crit][['disease', 'combo', 'gene', 'wrongFrac']]
        clinvar_final['panel_or_clinvar'] = 'Total ClinVar'
        clinvar_final.loc[:, 'gene'] = clinvar_final.apply(lambda row: row['gene'] + ' Total ClinVar', axis=1)

        df_final = pd.concat([clinvar_final, panel_final])
        df_final.loc[:, 'combo'] = df_final.apply(lambda row: row['combo'].replace('_base.', 'BASE ').replace('_trained.','TRAINED '), axis=1)
        df_final.to_csv(output.o, index=False, sep='\t')

rule paper_heatmap:
    input: WORK + 'paper_plot_data/heatmap'
    output: DOCS + 'paper_plts/fig5_heatmap.pdf'
    run:
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          p = ggplot(data=d) +
              geom_raster(aes(y=gene, x=reorder(combo, wrongFrac, mean), fill=wrongFrac)) +
              ylab('') + xlab('') + theme_bw(base_size=18) +
              scale_fill_gradient(low="yellow", high="blue") +
              labs(fill="Incorrect prediction fraction") +
              facet_grid(disease~., scale="free") +
              theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
          ggsave("{output}", p, width=15, height=9)
          """)

FIGS = ('fig1_count_plot', 'fig4_eval_clinvar', 'fig3_panelEval.byVarClassFalse', 'fig5_heatmap')
rule all_paper_plots:
    input: expand(DOCS + 'paper_plts/{fig}.pdf', fig=FIGS)
