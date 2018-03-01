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
