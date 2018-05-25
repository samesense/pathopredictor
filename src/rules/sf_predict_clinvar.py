"""Trained with panels,
   evaluate on clinvar
"""

def eval_pr(df, disease, acc_ls):
    scores = [x for x in df.columns.values if '_probaPred' in x or x in FEATS or x == 'revel']
    feat_names = {'Combination':'PathoPredictor', 'ccr':'CCR', 'fathmm':'FATHMM', 'revel':'REVEL',
                  'vest':'VEST', 'missense_badness':'Missense badness', 'missense_depletion':'Missense depletion'}
    for score in scores:
        fpr, tpr, _ = metrics.roc_curve(df['y'], df[score], pos_label=1)
        auc = metrics.auc(fpr, tpr)
        avg_pr = metrics.average_precision_score(df['y'], df[score])
        feature = score
        if '-' in feature:
            feature = 'Combination'
        acc_ls.append((auc, avg_pr, feat_names[feature], disease))

rule eval_by_gene_clinvar:
    input:  i = WORK + '{eval_set}/roc_df_clinvar/{features}'
    output: o = DATA + 'interim/EVAL_{eval_set}/pred_clinvar_eval/{eval_source}.{features}'
    run:
        df_pre = pd.read_csv(input.i, sep='\t')
        if wildcards.eval_source == 'clinvar_single':
            df = df_pre[df_pre.is_single]
        else:
            df = df_pre

        acc_ls = []
        for dis in set(df['Disease']):
            eval_pr(df[df.Disease==dis], dis, acc_ls)
        score_df = pd.DataFrame(acc_ls, columns=['auc', 'avg_pr', 'features', 'Disease'])
        benign_size_df = df[df.y==0].groupby(['Disease']).size().reset_index().rename(columns={0:'benign_size'})
        path_size_df = df[df.y==1].groupby(['Disease']).size().reset_index().rename(columns={0:'pathogenic_size'})
        size_df = pd.merge(benign_size_df, path_size_df, on='Disease', how='outer')
        m = pd.merge(score_df, size_df, on=['Disease'])
        clinvar_names  = {'clinvar_tot': 'Total ClinVar',
                          'clinvar_single': 'ClinVar w/ Evidence'}

        m.loc[:, 'clinvar_type'] = clinvar_names[wildcards.eval_source]
        disease_order = {'genedx-epi-limitGene':3,
                         'Rasopathies':4,
                         'EPI':2,
                         'Cardiomyopathy':1}

        diseases = {'genedx-epi-limitGene':'Epilepsy (dominant)',
                    'Rasopathies':'Rasopathies',
                    'EPI':'Epilepsy',
                    'Cardiomyopathy':'Cardiomyopathy'}
        m.loc[:, 'disease_name'] = m.apply(lambda row: diseases[row['Disease']], axis=1)
        m.loc[:, 'disease_order'] = m.apply(lambda row: disease_order[row['Disease']], axis=1)
        m.to_csv(output.o, index=False, sep='\t')

rule plot_clinvar_eval_paper:
    input:  expand(DATA + 'interim/EVAL_clinvar/pred_clinvar_eval/{clinvar_set}.' + C_FEATS, clinvar_set=('clinvar_tot', 'clinvar_single'))
    output: o = DOCS + 'paper_plts/fig6_evalClinvar.tiff'
    run:
        plot_cmd = """geom_col( aes(y=avg_pr, x=reorder(features, avg_pr)) ) +
                      geom_text(size=2, hjust="left", colour="white", data=label_df, aes(x=x, y=y, label=label))"""

        df_tot = pd.concat([pd.read_csv(afile, sep='\t') for afile in input])
        crit = df_tot.apply(lambda row: row['features'] != 'REVEL', axis=1)
        # df_tot.loc[:, 'feature_color'] = df_tot.apply(lambda row: 'bomdo' if row['features']=='Combination' else 'feat', axis=1)
        df_tot[crit].to_csv(output.o + '.df', index=False, sep='\t')

        df = df_tot[['clinvar_type', 'disease_name', 'benign_size', 'pathogenic_size']].drop_duplicates().melt(id_vars=['clinvar_type', 'disease_name'], var_name='var_type')
        df.loc[:, 'label'] = df.apply(lambda row: row['var_type'].split('_')[0][0] + '=%d' % (row['value']), axis=1)
        df.loc[:, 'x'] = df.apply(lambda row: 'Missense badness' if 'b' in row['label'] else 'Missense depletion', axis=1)
        df['y'] = .01
        df.to_csv(output.o + '.tmp.clinvar.labels', index=False, sep='\t')

        R("""
          require(ggplot2)
          feature_palette <- c("#D4ED91", "grey")
          d = read.delim("{output}.df", sep='\t', header=TRUE)
          label_df = read.delim("{output}.tmp.clinvar.labels", sep="\t", header=TRUE)
          d$clinvar_type = factor(d$clinvar_type, levels=c("Total ClinVar", "ClinVar w/ Evidence"))
          p = ggplot(data=d) + {plot_cmd} + guides(fill=FALSE) +
              ylab('Average precision') + xlab('') + theme_bw(base_size=10.5) + facet_grid(clinvar_type~disease_name) +
              coord_flip() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=10))
          ggsave("{output}", p, height=8, width=18, units="cm", dpi=300)
          """)
        shell('rm {output}.tmp.clinvar.labels')

rule plot_ndenovo_eval_paper:
    input:  i = DATA + 'interim/EVAL_ndenovo/pred_clinvar_eval/clinvar_tot.' + C_FEATS
    output: o = DOCS + 'paper_plts/fig8_evalDenovo.tiff'
    run:
        plot_cmd = """geom_col( aes(y=avg_pr, x=reorder(features, avg_pr)) )"""

        df_tot = pd.read_csv(input.i, sep='\t')
        crit = df_tot.apply(lambda row: 'Epi' in row['disease_name'] and ('REVEL'==row['features'] or 'PathoPredictor'==row['features']), axis=1)
        # df_tot.loc[:, 'feature_color'] = df_tot.apply(lambda row: 'bomdo' if row['features']=='Combination' else 'feat', axis=1)
        df_tot[crit].to_csv(output.o + '.df', index=False, sep='\t')

        R("""
          require(ggplot2)
          feature_palette <- c("#D4ED91", "grey")
          d = read.delim("{output}.df", sep='\t', header=TRUE)
          p = ggplot(data=d) + {plot_cmd} + guides(fill=FALSE) +
              ylab('Average precision') + xlab('') + theme_bw(base_size=12) + facet_grid(.~disease_name) +
              coord_flip() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=12))
          ggsave("{output}", p, height=4, width=12, units="cm", dpi=300)
          """)

rule both_clinvar_eval:
    input: DOCS + 'paper_plts/fig5_evalClinvar.tiff', DOCS + 'paper_plts/fig7_evalDenovo.tiff',


