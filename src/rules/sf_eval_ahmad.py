"""Combine evaluations by disease."""
import pandas as pd

#include: "const.py"
#include: "sf_eval_non_panel.py"

rule ahmad_percent_wrong:
    input:  panel = WORK + '{method}.eval_panel.{cols}.eval'
    output: o = WORK + '{method}.{cols}.eval_panel.eval.percentWrong'
    run:
        def calc_wrong_percent(rows):
            tot_wrong = list(rows[rows.eval_type == 'TotWrong']['var_count'].values)[0]
            tot_preds = sum(rows[rows.eval_type != 'TotWrong']['var_count'])
            return tot_wrong/tot_preds

        def calc_wrong_percent_benign(rows):
            tot_preds = sum(rows[ (rows.eval_type=='WrongBenign') | (rows.eval_type=='CorrectBenign')]['var_count'].values)
            ls = list(rows[rows.eval_type == 'WrongBenign']['var_count'].values)
            if ls:
                wrong_preds = ls[0]
            else:
                wrong_preds = 0
            return wrong_preds/tot_preds

        def calc_wrong_percent_path(rows):
            tot_preds = sum( rows[ (rows.eval_type=='WrongPath') | (rows.eval_type=='CorrectPath')]['var_count'].values )
            ls = list(rows[rows.eval_type == 'WrongPath']['var_count'].values)

            if ls:
                wrong_preds = ls[0]
            else:
                wrong_preds = 0

            return wrong_preds/tot_preds

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
#                if not 'no_disease' in score_type:
                return 'predict_' + score_type.split('::')[0] + '_' + wildcards.cols
                # else:
                #     prefix = 'panel_indep_'
                #     if '>2' in score_type:
                #         return prefix + 'paper_' + wildcards.cols
                #     return prefix + 'self-trained'               
                
            elif score_type == 'global_PredictionStatusMPC':
                return 'panel-trained_' + wildcards.cols
            elif 'global' == score_type[:6]:
                return score_type.replace('global_PredictionStatusMPC', 'train_with')
            return score_type

        def rename_disease(disease, score_type):
            if disease == 'no_disease':
                if '_holdOut' in score_type:
                    return 'self-trained'
                elif '>2' in score_type:
                    return 'paper-cutoff'
                else:
                    i = 1/0
            return disease

        # panel_df = (pd.read_csv(input.panel, sep='\t')
        #             .groupby(('disease', 'score_type'))
        #             .apply(calc_wrong_percent)
        #             .reset_index()
        #             .rename(columns={0:'percent_wrong'}) )

        panel_df_benign = (pd.read_csv(input.panel, sep='\t')
                           .groupby(('disease', 'score_type'))
                           .apply(calc_wrong_percent_benign)
                           .reset_index()
                           .rename(columns={0:'percent_wrong'}) )
        panel_df_benign['var_class'] = 'benign'
        panel_df_path= (pd.read_csv(input.panel, sep='\t')
                        .groupby(('disease', 'score_type'))
                        .apply(calc_wrong_percent_path)
                        .reset_index()
                        .rename(columns={0:'percent_wrong'}) )
        panel_df_path['var_class'] = 'pathogenic'

        panel_df = pd.concat([panel_df_path, panel_df_benign])
        panel_df['color'] = 'predict_panel'
        m = panel_df
        m.loc[:, 'st'] = m.apply(lambda row: rename_score(row['score_type']), axis=1)
        m.loc[:, 'dis'] = m.apply(lambda row: rename_disease(row['disease'], row['score_type']), axis=1)
        crit = m.apply(lambda row: not row['disease'] in ('ALL',)
                       and not 'uc' in row['disease']
                       and not 'earing' in row['disease']
                       and not 'issue' in row['disease']
                       and not ('no_disease' == row['disease'] and '>2' in row['score_type'])
                       and row['color'] != 'predict clinvar'
                       and not 'clinvar' in row['score_type'],
                       axis=1)
        m[crit].to_csv(output.o, index=False, sep='\t')

rule concat_extra:
    input:  expand( WORK + 'global.{cols}.eval_panel.eval.percentWrong', cols=COMBO_FEATS)
    output: o = WORK + 'cc'
    run:
        m = pd.concat([pd.read_csv(afile, sep='\t') for afile in list(input)]).drop_duplicates(subset=['color', 'disease', 'score_type', 'st', 'dis', 'var_class'])
        #print( m[['dis', 'percent_wrong']].groupby('dis').apply(min) )
        min_df = m[['dis', 'percent_wrong']].groupby('dis').apply(min).rename(columns={'dis':'dis_junk', 'percent_wrong':'min'}).reset_index()
        d = pd.merge(m, min_df, on='dis', how='left')
        d.loc[:, 'min_color'] = d.apply(lambda row: row['percent_wrong']==row['min'], axis=1)
        d.to_csv(output.o, sep='\t', index=False)

rule concat_extra_gene:
    input:  expand( WORK + 'global.{cols}.eval_panel.eval.percentWrong', cols=COMBO_FEATS)
    output: o = WORK + 'cc.gene'
    run:
        m = pd.concat([pd.read_csv(afile, sep='\t') for afile in list(input)]).drop_duplicates(subset=['color', 'disease', 'score_type', 'st', 'dis', 'var_class'])
        #print( m[['dis', 'percent_wrong']].groupby('dis').apply(min) )
        min_df = m[['dis', 'percent_wrong']].groupby('dis').apply(min).rename(columns={'dis':'dis_junk', 'percent_wrong':'min'}).reset_index()
        d = pd.merge(m, min_df, on='dis', how='left')
        d.loc[:, 'min_color'] = d.apply(lambda row: row['percent_wrong']==row['min'], axis=1)
        d.to_csv(output.o, sep='\t', index=False)
#theme(legend.position="none") 
rule plot_ahmad:
    input:  WORK + 'cc'
    output: DOCS + 'plot/global.byDisease.byVarClass{byVarClass}.png'
    run:
        if wildcards.byVarClass == 'True':
            plot_cmd = 'geom_col(aes(y=percent_wrong, x=reorder(st, percent_wrong), fill=var_class), position="dodge")'
        else:
            plot_cmd = 'geom_col(aes(y=percent_wrong, x=reorder(st, percent_wrong)), position="dodge")'
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          p = ggplot(data=d) + {plot_cmd} +
              facet_grid(dis~.) + theme_bw() +
              theme(axis.text.x = element_text(angle=90, hjust=1, size=8)) +
              ylab('Wrong prediction fraction') +
              xlab('') + coord_flip() + theme(axis.text.y = element_text(size=6))
          ggsave("{output}", p, height=20)
          """)

rule ahmad_prediction_plots:
    input: expand( DOCS + 'plot/global.byDisease.byVarClass{byVarClass}.png', byVarClass=('True', 'False') )
