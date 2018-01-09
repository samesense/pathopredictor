"""Combine evaluations by disease."""
import pandas as pd

include: "const.py"
include: "sf_eval_non_panel.py"

rule percent_wrong:
    input:  panel = WORK + '{method}.eval_panel.{cols}.eval',
            non_panel = expand( WORK + '{{method}}.eval_{dat}.{{cols}}.eval', \
                                dat=('clinvar', 'denovo', 'clinvar_mult', 'clinvar_single', 'clinvar_exp') )
    output: o = WORK + '{method}.{cols}.eval_panel.eval.percentWrong'
    run:
        def calc_wrong_percent(rows):
            print( len(rows) )
            tot_wrong = list(rows[rows.eval_type == 'TotWrong']['var_count'].values)[0]
            tot_preds = sum(rows[rows.eval_type != 'TotWrong']['var_count'])
            print('tw', tot_wrong)
            print('s', tot_preds)
            return tot_wrong/tot_preds

        def load_other_df(afile):
            panel_df = (pd.read_csv(afile, sep='\t')
                        .groupby(('clinvar_type', 'disease', 'score_type'))
                        .apply(calc_wrong_percent)
                        .reset_index()
                        .rename(columns={0:'percent_wrong', 'score_type':'score_type_pre'}) )
            panel_df['score_type'] = panel_df['clinvar_type'] + '::' + panel_df['score_type_pre']
#            print(panel_df.head())
            return panel_df[['disease', 'percent_wrong', 'score_type']]

        non_panel_dfs = pd.concat( [load_other_df(afile) for afile in input.non_panel] )
        panel_df = (pd.read_csv(input.panel, sep='\t')
                    .groupby(('disease', 'score_type'))
                    .apply(calc_wrong_percent)
                    .reset_index()
                    .rename(columns={0:'percent_wrong'})
                    .to_csv(output.o, index=False, sep='\t') )
        pd.concat([panel_df, non_panel_dfs]).to_csv(output.o, index=False, sep='\t')

rule plot_ahmad:
    input:  WORK + '{method}.{cols}.eval_panel.eval.percentWrong'
    output: DOCS + 'plot/{method}.{cols}.byDisease.png'
    run:
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          p = ggplot(data=d) +
          geom_col(aes(y=percent_wrong,x=score_type, fill=score_type)) +
          facet_grid(disease~., scale='free') + theme_bw() +
          ylab('Wrong Predictions') +
          theme(axis.text.x = element_text(angle=45, hjust=1)) +
          xlab('') + theme(legend.position="none")
          ggsave("{output}", p)
          """)

rule all_ahmad:
    input: expand( DOCS + 'plot/{method}.{cols}.byDisease.png', method=('global',), cols=('mpc', 'revel', 'mpc-revel') )
                   
    
