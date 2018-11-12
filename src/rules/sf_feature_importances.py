"""Compare feature importances for disease and panel vs clinvar"""

rule prep_varimpact:
    input:  i = DATA + 'interim/full/panel.dat'
    output: x = DATA + 'interim/varimpact_pre/panel.{disease}.x',
            y = DATA + 'interim/varimpact_pre/panel.{disease}.y'
    run:
        df = pd.read_csv(input.i, delimiter='\t')
        df.loc[:, 'ccr'] = df.apply(lambda row: row['ccr'] if row['ccr'] != 0 else random.uniform(0, 0.5), axis=1)
        df[df.Disease==wildcards.disease][FEATS + ['is_domain']].to_csv(output.x, sep='\t', header=True, index=False)
        df[df.Disease==wildcards.disease][['y']].to_csv(output.y, sep='\t', header=True, index=False)

rule varimpact:
    input:  DATA + 'interim/varimpact_pre/panel.{disease}.x',
            DATA + 'interim/varimpact_pre/panel.{disease}.y'
    output: DATA + 'interim/varimpact/panel.{disease}'
    singularity: "docker://samesense/varimpact-docker"
    shell:  'Rscript {SCRIPTS}varimpact.R {input} {output}'

rule tt:
    input: expand(DATA + 'interim/varimpact/panel.{d}', d=('EPI', 'Cardiomyopathy', 'Rasopathies',) )

rule feature_importance:
    input:  expand(DATA + 'interim/full/{eval_set}.dat', eval_set=('panel', 'clinvar') )
    output: DATA + 'interim/plot_data/importances'
    shell:  'python {SCRIPTS}feature_importance.py {input} {output}'

rule plot_feature_importance:
    input:  DATA + 'interim/plot_data/importances'
    output: DOCS + 'paper_plts/fig2_featureImportance.tiff'
    run:
        R("""
          require(ggplot2)
          d = read.delim("{input}", header=TRUE, sep="\t")
          dp = d[d$eval_set=="panel",]
          p = ggplot(data=dp, aes(y=importance, x=reorder(feature, importance, mean))) + \
              geom_col(colour="black", fill="gray") + \
              geom_errorbar(aes(ymin=importance-eb, ymax=importance+eb), size=.3, width=.2) + \
              facet_grid(disease~.) + coord_flip() + theme_bw(base_size=18) + \
              xlab('') + ylab('Feature importance')
          ggsave("{output}", p, dpi=300, height=22, width=10, units="cm")
          """)
# rule collapse_feature_importance:
#     input: expand(DATA + 'interim/importances/{eval_set}.feat_importance', eval_set=('panel',))
#     output: DATA + 'interim/plot_data/importance'
