"""Combine evaluations by disease."""

rule ahmad_percent_wrong:
    input:  panel = WORK + '{method}.eval_panel.{cols}.eval'
    output: o = WORK + '{method}.{cols}.eval_panel.eval.percentWrong'
    run:
        def calc_tot_vars(rows):
            tot_preds = sum(rows[rows.eval_type != 'TotWrong']['var_count'])
            return tot_preds

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

        panel_size_df = (pd.read_csv(input.panel, sep='\t')
                         .groupby(('disease', 'score_type'))
                         .apply(calc_tot_vars)
                         .reset_index()
                         .rename(columns={0:'tot_vars'}) )

        panel_wrong_df = (pd.read_csv(input.panel, sep='\t')
                          .groupby(('disease', 'score_type'))
                          .apply(calc_wrong_percent)
                          .reset_index()
                          .rename(columns={0:'percent_wrong'}) )
        panel_df = pd.merge(panel_wrong_df, panel_size_df, on=['disease', 'score_type'])

        # panel_df_benign = (pd.read_csv(input.panel, sep='\t')
        #                    .groupby(('disease', 'score_type'))
        #                    .apply(calc_wrong_percent_benign)
        #                    .reset_index()
        #                    .rename(columns={0:'percent_wrong'}) )
        # panel_df_benign['var_class'] = 'benign'
        # panel_df_path= (pd.read_csv(input.panel, sep='\t')
        #                 .groupby(('disease','score_type'))
        #                 .apply(calc_wrong_percent_path)
        #                 .reset_index()
        #                 .rename(columns={0:'percent_wrong'}) )
        # panel_df_path['var_class'] = 'pathogenic'

        #panel_df = pd.concat([panel_df_path, panel_df_benign])
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

def color_bar(row):
    if row['st'] in ('paper_mpc', 'paper_revel', 'paper_ccr'):
        return 'Baseline'
    if 'paper' in row['st']:
        return 'Combined baseline'
    return 'Trained'

rule concat_extra:
    input:  expand( WORK + 'global.{cols}.eval_panel.eval.percentWrong', cols=COMBO_FEATS)
    output: o = WORK + 'cc'
    run:
        m = pd.concat([pd.read_csv(afile, sep='\t') for afile in list(input)]).drop_duplicates(subset=['color', 'disease', 'score_type', 'st', 'dis'])
        #print( m[['dis', 'percent_wrong']].groupby('dis').apply(min) )
        min_df = m[['dis', 'percent_wrong']].groupby('dis').apply(min).rename(columns={'dis':'dis_junk', 'percent_wrong':'min'}).reset_index()
        dp = pd.merge(m, min_df, on='dis', how='left')
        disease_order = {'genedx-epi-limitGene':3,
                         'Rasopathies':4,
                         'genedx-epi':2,
                         'Cardiomyopathy':1}

        diseases = {'genedx-epi-limitGene':'Epilepsy (dominant genes)',
                    'Rasopathies':'Rasopathies',
                    'genedx-epi':'Epilepsy',
                    'Cardiomyopathy':'Cardiomyopathy'}
        crit = dp.apply(lambda row: row['dis'] in diseases, axis=1)
        d = dp[crit]
        d.loc[:, 'dis_order'] = d.apply(lambda row: disease_order[row['dis']], axis=1)
        d.loc[:, 'dis'] = d.apply(lambda row: diseases[row['dis']] + ' {n=%d}'
                                  % (row['tot_vars']), axis=1)
        d.loc[:, 'Classifier'] = d.apply(color_bar, axis=1)
        d.loc[:, 'is_best'] = d.apply(lambda row: row['percent_wrong']==row['min'], axis=1)
        d.loc[:, 'st'] = d.apply(lambda row: row['st'].replace('panel-trained_','TRAINED_').replace('paper_','BASE_'), axis=1)
        d.to_csv(output.o, sep='\t', index=False)

rule concat_extra_gene:
    input:  expand( WORK + 'global.{cols}.eval_panel.eval.percentWrong', cols=COMBO_FEATS)
    output: o = WORK + 'cc.gene'
    run:
        m = pd.concat([pd.read_csv(afile, sep='\t') for afile in list(input)]).drop_duplicates(subset=['color', 'disease', 'score_type', 'st', 'dis', 'var_class'])
        min_df = m[['dis', 'percent_wrong']].groupby('dis').apply(min).rename(columns={'dis':'dis_junk', 'percent_wrong':'min'}).reset_index()
        d = pd.merge(m, min_df, on='dis', how='left')
        d.loc[:, 'min_color'] = d.apply(lambda row: row['percent_wrong']==row['min'], axis=1)
        d.loc[:, 'st'] = d.apply(lambda row: row['st'].replace('panel-trained_','').replace('paper_',''), axis=1)
        d.to_csv(output.o, sep='\t', index=False)

rule plot_ahmad:
    input:  WORK + 'cc'
    output: DOCS + 'paper_plts/fig3_panelEval.byVarClass{byVarClass}.pdf'
    run:
        if wildcards.byVarClass == 'True':
            plot_cmd = 'geom_col(aes(y=percent_wrong, x=reorder(st, percent_wrong), colour=is_best, fill=classifier_color), position="dodge")'
        else:
            plot_cmd = 'geom_col(aes(fill=Classifier, y=percent_wrong, x=reorder(st, percent_wrong))) + geom_point(data=dbest, aes(x=st,y=percent_wrong))'
        R("""
          require(ggplot2)
          d = read.delim("{input}", sep='\t', header=TRUE)
          d$dis = factor(d$dis, levels=unique( d[order(d$dis_order),]$dis ))
          dbest = d[d$is_best=="True",]
          head(dbest)
          p = ggplot(data=d) + {plot_cmd} +
              facet_grid(.~dis) + theme_bw() +
              theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=12)) +
              ylab('Incorrect prediction fraction') + theme(legend.position="bottom") +
              xlab('') + coord_flip() + theme(axis.text.y = element_text(size=10))
          ggsave("{output}", p, width=20)
          """)

rule ahmad_prediction_plots:
    input: expand( DOCS + 'paper_plts/fig3_panelEval.byVarClass{byVarClass}.pdf', byVarClass=('False',) )
