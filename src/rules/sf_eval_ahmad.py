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
            tot_wrong = list(rows[rows.eval_type == 'TotWrong']['var_count'].values)[0]
            tot_preds = sum(rows[rows.eval_type != 'TotWrong']['var_count'])
            # print('tw', tot_wrong)
            # print('s', tot_preds)
            return tot_wrong/tot_preds

        def load_other_df(afile):
            panel_df = (pd.read_csv(afile, sep='\t')
                        .groupby(('clinvar_type', 'disease', 'score_type'))
                        .apply(calc_wrong_percent)
                        .reset_index()
                        .rename(columns={0:'percent_wrong', 'score_type':'score_type_pre'}) )
            panel_df['st'] = panel_df.apply(lambda row: row['score_type_pre'].split(row['disease'])[0], axis=1)
            panel_df['score_type'] = panel_df['clinvar_type'] + '::' + panel_df['st']
            return panel_df[['disease', 'percent_wrong', 'score_type']]

        def rename_score(score_type):
            if score_type == 'global_PredictionStatusMPC>2':
                return 'paper_' + wildcards.cols
            elif '::' in score_type:
                return 'predict_' + score_type.split('::')[0]
            elif score_type == 'global_PredictionStatusMPC':
                return 'panel-trained'
            elif 'global' == score_type[:6]:
                return score_type.replace('global_PredictionStatusMPC', 'train_with')
            return score_type            
                
        non_panel_dfs = pd.concat( [load_other_df(afile) for afile in input.non_panel] )
        non_panel_dfs['color'] = 'predict clinvar'
        panel_df = (pd.read_csv(input.panel, sep='\t')
                    .groupby(('disease', 'score_type'))
                    .apply(calc_wrong_percent)
                    .reset_index()
                    .rename(columns={0:'percent_wrong'}) )
        panel_df['color'] = 'predict_panel'
        m = pd.concat([panel_df, non_panel_dfs])
        m.loc[:, 'st'] = m.apply(lambda row: rename_score(row['score_type']), axis=1)
        crit = m.apply(lambda row: not row['disease'] in ('ALL', 'no_disease')
                       and not 'earing' in row['disease']
                       and not 'issue' in row['disease'],
                       axis=1)
        m[crit].to_csv(output.o, index=False, sep='\t')

rule plot_ahmad:
    input:  WORK + '{method}.{cols}.eval_panel.eval.percentWrong'
    output: DOCS + 'plot/{method}.{cols}.byDisease.png'
    run:
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          p = ggplot(data=d) +
          geom_col(aes(y=percent_wrong, x=st, fill=color)) +
          facet_grid(disease~.) + theme_bw() +
          theme(axis.text.x = element_text(angle=90, hjust=1)) +
          ylab('Wrong prediction fraction') +
          xlab('') + theme(legend.position="none") + coord_flip()
          ggsave("{output}", p, height=10)
          """)

rule all_ahmad:
    input: expand( DOCS + 'plot/{method}.{cols}.byDisease.png', method=('global',), cols=('mpc', 'revel', 'mpc-revel') )
                   
    
