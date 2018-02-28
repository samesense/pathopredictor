"""Predictor paper plots."""

rule count_plot_data:
    input: panel = WORK + 'roc_df_panel/mpc',
           clinvar = WORK + 'roc_df_clinvar/mpc'
    output: o = WORK + 'paper_plot_data/count_plot'
    run:
        diseases = set(('Rasopathies', 'genedx-epi', 'Cardiomyopathy'))
        df_panel = pd.read_csv(input.panel, sep='\t')
        crit = df_panel.apply(lambda row: row['Disease'] in diseases, axis=1)
        panel = df_panel[crit].groupby(['Disease','y']).size().reset_index().rename(columns={0:'var_count'})
        panel['eval_type'] = 'panel'
        df_clinvar = pd.read_csv(input.clinvar, sep='\t')
        crit = df_clinvar.apply(lambda row: row['Disease'].split(':')[0] in diseases and row['Disease'].split(':')[1] in ('tot', 'single'), axis=1)
        clinvar = df_clinvar[crit].groupby(['Disease','y']).size().reset_index().rename(columns={0:'var_count'})
        clinvar.loc[:, 'eval_type'] = clinvar.apply(lambda row: 'ClinVar ' + row['Disease'].split(':')[1], axis=1)
        clinvar.loc[:, 'Disease'] = clinvar.apply(lambda row: row['Disease'].split(':')[0], axis=1)
        df = pd.concat([panel, clinvar])
        df.to_csv(output.o, index=False, sep='\t')

rule count_plot:
    input:  WORK + 'paper_plot_data/count_plot'
    output: DOCS + 'paper_plts/fig1_count_plot.pdf'
    run:
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          p = ggplot(data=d) +
              geom_bar(stat="identity", aes(x=Disease,y=var_count,fill=factor(y))) +
              facet_grid(eval_type~.) + theme_bw(base_size=18)
          ggsave("{output}", p)
          """)
