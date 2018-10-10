"""Trained with panels,
   evaluate on clinvar
"""
def eval_pr_curve_clinvar(df, disease, acc_ls, out, clinvar_type, use_revel=True):
    """Dump pr curve data."""
    scores = [x for x in df.columns.values if '_probaPred' in x or x in FEATS or x == 'revel']
    if not use_revel:
        scores = [x for x in df.columns.values if '_probaPred' in x or x in FEATS]
    feat_names = {'Combination':'PathoPredictor', 'ccr':'CCR', 'fathmm':'FATHMM', 'revel':'REVEL',
                  'vest':'VEST', 'missense_badness':'Missense badness', 'missense_depletion':'Missense depletion'}
    disease_order = {'genedx-epi-limitGene':3,
                         'Rasopathies':4,
                         'EPI':2,
                         'Cardiomyopathy':1}

    diseases = {'genedx-epi-limitGene':'Epilepsy (dominant)',
                    'Rasopathies':'Rasopathies',
                    'EPI':'Epilepsy',
                    'Cardiomyopathy':'Cardiomyopathy'}

    for score in scores:
        precision, recall, thresholds = precision_recall_curve(df['y'], df[score], pos_label=1)
        feature = score
        if '-' in feature:
            feature = 'Combination'
        disease_name = diseases[disease]
        order = disease_order[disease]
        for pre, rec in zip(precision, recall):
            ls = [str(x) for x in (feat_names[feature], disease, pre, rec, disease_name, order, clinvar_type,) ]
            print('\t'.join(ls), file=out)


def eval_pr_curve(df, disease, acc_ls, out, use_revel=True):
    """Dump pr curve data."""
    scores = [x for x in df.columns.values if '_probaPred' in x or x in FEATS or x == 'revel']
    if not use_revel:
        scores = [x for x in df.columns.values if '_probaPred' in x or x in FEATS]
    feat_names = {'Combination':'PathoPredictor', 'ccr':'CCR', 'fathmm':'FATHMM', 'revel':'REVEL',
                  'vest':'VEST', 'missense_badness':'Missense badness', 'missense_depletion':'Missense depletion'}
    disease_order = {'genedx-epi-limitGene':3,
                         'Rasopathies':4,
                         'EPI':2,
                         'Cardiomyopathy':1}

    diseases = {'genedx-epi-limitGene':'Epilepsy (dominant)',
                    'Rasopathies':'Rasopathies',
                    'EPI':'Epilepsy',
                    'Cardiomyopathy':'Cardiomyopathy'}

    for score in scores:
        precision, recall, thresholds = precision_recall_curve(df['y'], df[score], pos_label=1)
        feature = score
        if '-' in feature:
            feature = 'Combination'
        disease_name = diseases[disease]
        order = disease_order[disease]
        for pre, rec in zip(precision, recall):
            ls = [str(x) for x in (feat_names[feature], disease, pre, rec, disease_name, order)]
            print('\t'.join(ls), file=out)


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
        if feature in feat_names:
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

rule plot_clinvar_eval_pr:
    input:  i = WORK + 'clinvar/roc_df_clinvar/{features}'
    output: o = DATA + 'interim/EVAL_{eval_set}/pred_clinvar_eval_curve/{eval_source}.{features}'
    run:
        df_pre = pd.read_csv(input.i, sep='\t')
        if wildcards.eval_source == 'clinvar_single':
            df = df_pre[df_pre.is_single]
        else:
            df = df_pre

        clinvar_names  = {'clinvar_tot': 'Total ClinVar',
                          'clinvar_single': 'ClinVar w/ Evidence'}


        acc_ls = []
        with open(output.o, 'w') as fout:
            print('features\tDisease\tPrecision\tRecall\tdisease_name\tdisease_order\tclinvar_type', file=fout)
            for dis in set(df['Disease']):
                eval_pr_curve_clinvar(df[df.Disease==dis], dis, acc_ls, fout, clinvar_names[wildcards.eval_source], use_revel=True)
 
rule plot_clinvar_pr_curve:
    input:  i = expand(DATA + 'interim/EVAL_clinvar/pred_clinvar_eval_curve/{clinvar_set}.' + C_FEATS, clinvar_set=('clinvar_tot', 'clinvar_single'))
    output: o = DOCS + 'paper_plts/fig6a_curve.tiff'
    run:
        plot_cmd = """geom_line( aes(y=Precision, x=Recall,colour=features, group=features, fill=features))"""
        df_tot = pd.concat([pd.read_csv(afile, sep='\t') for afile in input])
        df_tot.to_csv(output.o + '.df', index=False, sep='\t')

        R("""
          require(ggplot2)
          require(grid)
          d = read.delim("{output}.df", sep='\t', header=TRUE)
          d$disease_name = factor(d$disease_name, levels=unique( d[order(d$disease_order),]$disease_name ))
          p = ggplot(data=d) + {plot_cmd} +
              facet_grid(clinvar_type~disease_name) + theme_bw(base_size=10) +
              theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=10)) +
              ylab('Precision') + labs(colour = "", fill="") +
              xlab('Recall')
          tiff("{output}", res=300, units="cm", height=7.5, width=18)
          grid.draw(p)
          grid.text("a", x=0.05, y=0.96)
          dev.off()
          """)

rule plot_clinvar_eval_paper:
    input:  expand(DATA + 'interim/EVAL_clinvar/pred_clinvar_eval/{clinvar_set}.' + C_FEATS, clinvar_set=('clinvar_tot', 'clinvar_single'))
    output: o = DOCS + 'paper_plts/fig6b_bar.tiff'
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
          require(grid)
          feature_palette <- c("#D4ED91", "grey")
          d = read.delim("{output}.df", sep='\t', header=TRUE)
          label_df = read.delim("{output}.tmp.clinvar.labels", sep="\t", header=TRUE)
          d$clinvar_type = factor(d$clinvar_type, levels=c("Total ClinVar", "ClinVar w/ Evidence"))
          p = ggplot(data=d) + {plot_cmd} + guides(fill=FALSE) +
              ylab('Average precision') + xlab('') + theme_bw(base_size=10.5) + facet_grid(clinvar_type~disease_name) +
              coord_flip() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=10))
          tiff("{output}", res=300, units="cm", height=7.5, width=18)
          grid.draw(p)
          grid.text("b", x=0.05, y=0.96)
          dev.off()
          """)
        shell('rm {output}.tmp.clinvar.labels')

rule plot_ndenovo_pr_curve_paper:
    input:  i = DATA + 'interim/EVAL_ndenovo/pred_clinvar_eval_curve/clinvar_tot.' + C_FEATS
    output: o = DOCS + 'paper_plts/fig7b_evalDenovoCurve.tiff'
    run:
        plot_cmd = """geom_line( aes(y=Precision, x=Recall,colour=features, group=features, fill=features))"""
        df_tot = pd.read_csv(input.i, sep='\t')
        crit = df_tot.apply(lambda row: 'Epi' in row['disease_name'] and ('REVEL'==row['features'] or 'PathoPredictor'==row['features']), axis=1)
        df_tot[crit].to_csv(output.o + '.df', index=False, sep='\t')

        R("""
          require(ggplot2)
          d = read.delim("{output}.df", sep='\t', header=TRUE)
          d$disease_name = factor(d$disease_name, levels=unique( d[order(d$disease_order),]$disease_name ))
          p = ggplot(data=d) + {plot_cmd} +
              facet_grid(.~disease_name) + theme_bw(base_size=10) +
              theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=10)) +
              ylab('Precision') + labs(colour = "", fill="") +
              xlab('Recall')
          ggsave("{output}", p, dpi=300, width=10, height=3.5, units="cm")
          """)

        # R("""
        #   require(ggplot2)
        #   require(grid)
        #   feature_palette <- c("#D4ED91", "grey")
        #   d = read.delim("{output}.df", sep='\t', header=TRUE)
        #   p = ggplot(data=d) + {plot_cmd} + guides(fill=FALSE) +
        #       ylab('Average precision') + xlab('') + theme_bw(base_size=12) + facet_grid(.~disease_name) +
        #       coord_flip() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=12)) 
        #   tiff("{output}", res=300, units="cm", height=3.5, width=10)
        #   grid.draw(p)
        #   grid.text("b", x=0.05, y=0.96)
        #   dev.off()
        #   """)

rule plot_ndenovo_eval_paper:
    input:  i = DATA + 'interim/EVAL_ndenovo/pred_clinvar_eval/clinvar_tot.' + C_FEATS
    output: o = DOCS + 'paper_plts/fig7c_evalDenovo.tiff'
    run:
        plot_cmd = """geom_col( aes(y=avg_pr, x=reorder(features, avg_pr)) )"""

        df_tot = pd.read_csv(input.i, sep='\t')
        crit = df_tot.apply(lambda row: 'Epi' in row['disease_name'] and ('REVEL'==row['features'] or 'PathoPredictor'==row['features']), axis=1)
        # df_tot.loc[:, 'feature_color'] = df_tot.apply(lambda row: 'bomdo' if row['features']=='Combination' else 'feat', axis=1)
        df_tot[crit].to_csv(output.o + '.df', index=False, sep='\t')

        R("""
          require(ggplot2)
          require(grid)
          feature_palette <- c("#D4ED91", "grey")
          d = read.delim("{output}.df", sep='\t', header=TRUE)
          p = ggplot(data=d) + {plot_cmd} + guides(fill=FALSE) +
              ylab('Average precision') + xlab('') + theme_bw(base_size=12) + facet_grid(.~disease_name) +
              coord_flip() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=12)) 
          tiff("{output}", res=300, units="cm", height=3.5, width=10)
          grid.draw(p)
          grid.text("b", x=0.05, y=0.96)
          dev.off()
          """)

rule combine_figs_train_panel_test_clinvar:
    input:  DOCS + 'paper_plts/fig6a_curve.tiff',
            DOCS + 'paper_plts/fig6b_bar.tiff'
    output: o = DOCS + 'paper_plts/fig6_evalClinvar.tiff'
    singularity:
        'docker://ncsapolyglot/converters-imagemagick'
    shell:  'convert -append {input} {output}'

