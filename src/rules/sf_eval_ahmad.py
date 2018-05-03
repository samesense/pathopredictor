"""Combine evaluations by disease."""

rule mk_panel_eval_figure_data:
    input:  i = WORK + 'roc_df_panel/{features}'
    output: o = DATA + 'interim/pred_panel_eval/{features}'
    run:
        df = pd.read_csv(input.i, sep='\t')
        acc_ls = []
        for dis in set(df['Disease']):
            eval_pr(df[df.Disease==dis], dis, acc_ls)
        score_df = pd.DataFrame(acc_ls, columns=['auc', 'avg_pr', 'features', 'Disease'])
        benign_size_df = df[df.y==0].groupby(['Disease']).size().reset_index().rename(columns={0:'benign_size'})
        path_size_df = df[df.y==1].groupby(['Disease']).size().reset_index().rename(columns={0:'pathogenic_size'})
        size_df = pd.merge(benign_size_df, path_size_df, on='Disease', how='outer')
        m = pd.merge(score_df, size_df, on=['Disease'])

        disease_order = {'genedx-epi-limitGene':3,
                         'Rasopathies':4,
                         'EPI':2,
                         'Cardiomyopathy':1}

        diseases = {'genedx-epi-limitGene':'Epilepsy (dominant genes)',
                    'Rasopathies':'Rasopathies',
                    'EPI':'Epilepsy',
                    'Cardiomyopathy':'Cardiomyopathy'}
        m.loc[:, 'disease_name'] = m.apply(lambda row: diseases[row['Disease']], axis=1)
        m.loc[:, 'disease_order'] = m.apply(lambda row: disease_order[row['Disease']], axis=1)
        m.to_csv(output.o, index=False, sep='\t')

rule plot_panel_eval:
    input:  i = DATA + 'interim/pred_panel_eval/' + C_FEATS
    output: o = DOCS + 'paper_plts/fig3_panelEval.pdf'
    run:
        plot_cmd = """geom_col( aes(y=avg_pr, x=reorder(features, avg_pr)) ) +
                      geom_text(data=label_df, colour="white", aes(x=x, y=y, label=label), hjust="left")"""
        df = pd.read_csv(input.i, sep='\t')[['disease_name', 'benign_size', 'pathogenic_size']].drop_duplicates().melt(id_vars=['disease_name'], var_name='var_type')
        df.loc[:, 'label'] = df.apply(lambda row: row['var_type'].split('_')[0] + '=%d' % (row['value']), axis=1)
        df.loc[:, 'x'] = df.apply(lambda row: 'Missense badness' if 'benign' in row['label'] else 'FATHMM', axis=1)
        df['y'] = 0.01
        df.to_csv(output.o + '.tmp.panel.labels', index=False, sep='\t')

        R("""
          require(ggplot2)
          label_df = read.csv('{output}.tmp.panel.labels', sep='\t')
          d = read.delim("{input}", sep='\t', header=TRUE)
          d$disease_name = factor(d$disease_name, levels=unique( d[order(d$disease_order),]$disease_name ))
          p = ggplot(data=d) + {plot_cmd} +
              facet_grid(.~disease_name) + theme_bw(base_size=18) +
              theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=12)) +
              ylab('Average precision') +
              xlab('') + coord_flip() + theme(axis.text.y = element_text(size=10))
          ggsave("{output}", p, width=20)
          """)
        shell('rm {output}.tmp.panel.labels')
