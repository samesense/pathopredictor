"""Combine evaluations by disease."""

rule auc_roc_and_avg_pre_anova:
    input:  i = WORK + 'roc_df_{data}/{features}'
    output: o = WORK + 'eval_features_{data}/{features}'
    shell:  'python {SCRIPTS}score_features.py {wildcards.features} {input} {output}'

rule combine_auc_avgPre_improveProb:
    input:  auc = WORK + 'eval_features_{data}/{features}',
            ip = DATA + 'interim/improveProb_out_collapse/{features}'
    output: o = DATA + 'interim/fig3_data_{data}/{features}'
    run:
        pd.merge(pd.read_csv(input.auc, sep='\t'),
                 pd.read_csv(input.ip, sep='\t'),
                 on=['Disease', 'combo'], how='left').to_csv(output.o, index=False, sep='\t')

rule combine_improveProb_features:
    input: expand(DATA + 'interim/fig3_data_panel/{feats}', feats=COMBO_FEATS)
    output: o = DATA + 'interim/fig3_data_panel.improveProb'
    run:
        pd.concat([pd.read_csv(afile, sep='\t') for afile in input]).to_csv(output.o, index=False, sep='\t')

rule ahmad_percent_wrong:
    input:  panel = WORK + '{method}.eval_panel.{cols}.eval'
    output: o = WORK + '{method}.{cols}.eval_panel.eval.percentWrong'
    run:
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
                         .reset_index() )

        panel_wrong_path_df = (pd.read_csv(input.panel, sep='\t')
                               .groupby(('disease', 'score_type'))
                               .apply(calc_wrong_percent_path)
                               .reset_index()
                               .rename(columns={0:'percent_wrong_path'}) )
        panel_wrong_benign_df = (pd.read_csv(input.panel, sep='\t')
                                 .groupby(('disease', 'score_type'))
                                 .apply(calc_wrong_percent_benign)
                                 .reset_index()
                                 .rename(columns={0:'percent_wrong_benign'}) )
        panel_wrong_df = pd.merge(panel_wrong_path_df, panel_wrong_benign_df, on=('disease', 'score_type'), how='outer')
        panel_wrong_df['percent_wrong'] = (panel_wrong_df['percent_wrong_path'] +
                                           panel_wrong_df['percent_wrong_benign'])/2

        panel_df = pd.merge(panel_wrong_df, panel_size_df, on=['disease', 'score_type'])

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
    if row['st'] in ('paper_mpc', 'paper_revel', 'paper_ccr', 'paper_domain'):
        return 'Baseline'
    if 'paper' in row['st']:
        return 'Combined baseline'
    return 'Trained'

def read_it(afile):
    combo = afile.split('/')[-1].split('.')[1]
    df = pd.read_csv(afile, sep='\t')
    df['combo'] = combo
    return df

rule concat_extra:
    input:  expand( WORK + 'global.{cols}.eval_panel.eval.percentWrong', cols=COMBO_FEATS)
    output: o = WORK + 'cc'
    run:
        m = pd.concat([read_it(afile) for afile in input]).drop_duplicates(subset=['color', 'disease', 'score_type', 'st', 'dis'])
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
        d.loc[:, 'dis'] = d.apply(lambda row: diseases[row['dis']], axis=1)
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

rule panel_pvals:
    """Test each trained against baseline w/ fishers test"""
    input:  WORK + 'cc',
            DATA + 'interim/fig3_data_panel.improveProb'
    output: o = WORK + 'cc.pvals'
    shell:  'python {SCRIPTS}mk_panel_pvalues.py {input} {output}'

rule plot_ahmad:
    input:  i = WORK + 'cc.pvals'
    output: DOCS + 'paper_plts/fig3_panelEval.byVarClass{byVarClass}.pdf'
    run:
        if wildcards.byVarClass == 'True':
            plot_cmd = 'geom_col(aes(y=percent_wrong, x=reorder(st, percent_wrong), colour=is_best, fill=classifier_color), position="dodge")'
        else:
            plot_cmd = """geom_col(aes(fill=Classifier, y=percent_wrong, x=reorder(st, percent_wrong, median))) +
                          geom_text(data=label_df, aes(x=x1,y=y,label=label_path), hjust=0) +
                          geom_text(data=label_df, aes(x=x2,y=y,label=label_benign), hjust=0)"""
        df = pd.read_csv(input.i, sep='\t')[['dis','tot_vars','CorrectPath','WrongPath','CorrectBenign','WrongBenign']].drop_duplicates()
        df.loc[:, 'label_path'] = df.apply(lambda row: 'pathogenic=%d' % (row['CorrectPath'] + row['WrongPath']), axis=1)
        df.loc[:, 'label_benign'] = df.apply(lambda row: 'benign=%d' % (row['CorrectBenign'] + row['WrongBenign']), axis=1)
        df['y'] = 0.4
        df['x1'] = 'TRAINED_mpc-revel-ccr'
        df['x2'] = 'TRAINED_revel-ccr'
        df.to_csv('tmp.labels', index=False, sep='\t')
        R("""
          require(ggplot2)
          # The palette with grey: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
          cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
          box_colors = c('white', 'black')
          label_df = read.csv('tmp.labels', sep='\t')
          d = read.delim("{input}", sep='\t', header=TRUE)
          d$dis = factor(d$dis, levels=unique( d[order(d$dis_order),]$dis ))
          dbest = d[d$is_best=="True",]
          head(dbest)
          p = ggplot(data=d) + {plot_cmd} +
              facet_grid(.~dis) + theme_bw(base_size=18) +
              theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=12)) +
              ylab('Incorrect prediction fraction') + theme(legend.position="bottom") +
              xlab('') + coord_flip() + theme(axis.text.y = element_text(size=10)) + scale_fill_manual(values=cbPalette)
          ggsave("{output}", p, width=20)
          """)
        shell('rm tmp.labels')

rule plot_idi:
    input:  i = WORK + 'cc.pvals'
    output: o = DOCS + 'paper_plts/fig4_idi.pdf'
    run:
        df = pd.read_csv(input.i, sep='\t')[['dis','tot_vars','CorrectPath','WrongPath','CorrectBenign','WrongBenign']].drop_duplicates()
        df.loc[:, 'label_path'] = df.apply(lambda row: 'pathogenic=%d' % (row['CorrectPath'] + row['WrongPath']), axis=1)
        df.loc[:, 'label_benign'] = df.apply(lambda row: 'benign=%d' % (row['CorrectBenign'] + row['WrongBenign']), axis=1)
        df['y'] = 0.2
        df['x1'] = 'TRAINED_mpc-ccr'
        df['x2'] = 'TRAINED_ccr-is_domain'
        tmp_labels = output.o + 'tmplabels'
        df.to_csv(tmp_labels, index=False, sep='\t')

        df = pd.read_csv(input.i, sep='\t')
        crit = df.apply(lambda row: not 'BASE' in row['st'] and '-' in row['st'], axis=1)
        tmp = output.o + '.tmp'
        df[crit].to_csv(tmp, index=False, sep='\t')

        plot_cmd = """geom_col(fill="#56B4E9", aes(alpha=box, y=idi, x=reorder(st, idi, median))) + scale_alpha_manual(values=c(.5, 1), guide="none") + geom_errorbar(aes(x=reorder(st, idi, median), ymin=idi_lower, ymax=idi_upper), width=.2, position=position_dodge(.9))"""
        #plot_cmd = """geom_col(fill="#56B4E9", aes(colour=box, y=idi, x=reorder(st, idi, median)))"""

        R("""
          require(ggplot2)
          cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
          box_colors = c('white', 'black')
          label_df = read.csv('{tmp_labels}', sep='\t')
          d = read.delim("{tmp}", sep='\t', header=TRUE)
          d$dis = factor(d$dis, levels=unique( d[order(d$dis_order),]$dis ))
          p = ggplot(data=d) + {plot_cmd} +
              facet_grid(.~dis) + theme_bw(base_size=18) +
              theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=14)) +
              ylab('Integrated discrimination index') + theme(legend.position="none") +
              xlab('') + coord_flip() + theme(axis.text.y = element_text(size=12)) + scale_colour_manual(values=box_colors)
          ggsave("{output}", p, width=20)

        """)
        shell('rm {tmp_labels}')
        shell('rm {tmp}')
        
rule ahmad_prediction_plots:
    input: expand( DOCS + 'paper_plts/fig3_panelEval.byVarClass{byVarClass}.pdf', byVarClass=('False',) )
