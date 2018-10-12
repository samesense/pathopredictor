rule plot_ndenovo_eval_paper:
    input:  i = DATA + 'interim/EVAL_ndenovo/pred_clinvar_eval/clinvar_tot.' + C_FEATS
    output: o = DOCS + 'paper_plts/fig8c_evalDenovoAvgPr.tiff'
    run:
        plot_cmd = """geom_col( aes(y=avg_pr, x=reorder(features, avg_pr)) )"""

        df_tot = pd.read_csv(input.i, sep='\t')
        crit = df_tot.apply(lambda row: 'Epi' in row['disease_name'] and ('MPC'==row['features'] or 'REVEL'==row['features'] or 'PathoPredictor'==row['features']), axis=1)
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
          grid.text("c", x=0.05, y=0.96)
          dev.off()
          """)

rule plot_ndenovo_pr_curve_paper:
    input:  i = DATA + 'interim/EVAL_ndenovo/pred_clinvar_eval_curve/clinvar_tot.' + C_FEATS
    output: o = DOCS + 'paper_plts/fig8b_evalDenovoCurve.tiff'
    run:
        plot_cmd = """geom_line( aes(y=Precision, x=Recall, colour=new_feat_name))"""
        df_tot = pd.read_csv(input.i, sep='\t')
        #crit = df_tot.apply(lambda row: 'Epi' in row['disease_name'] and 'MPC'==row['features'], axis=1)
        #crit = df_tot.apply(lambda row: 'Epi' in row['disease_name'] and 'PathoPredictor'==row['features'], axis=1)
        crit = df_tot.apply(lambda row: 'Epi' in row['disease_name'] and ('MPC' == row['features'] or 'REVEL'==row['features'] or 'PathoPredictor'==row['features']), axis=1)
        #crit = df_tot.apply(lambda row: 'Epi' in row['disease_name'] and ('MPC' == row['features'] or 'PathoPredictor'==row['features']), axis=1)
        #crit = df_tot.apply(lambda row: 'Epi' in row['disease_name'] and ('REVEL'==row['features'] or 'PathoPredictor'==row['features']), axis=1)
        plot_df = df_tot[crit]
        plot_df.loc[:, 'new_feat_name'] = plot_df.apply(lambda row: row['features'] if 'PathoPredictor' != row['features'] else row['features'] + '\n' + row['disease_name'], axis=1)
        crit = plot_df.apply(lambda row: 'Epilepsy' == row['disease_name'] or 'PathoPredictor'==row['features'], axis=1)
        plot_df[crit].to_csv(output.o + '.df', index=False, sep='\t')
        R("""
          require(ggplot2)
          require(grid)
          d = read.delim("{output}.df", sep='\t', header=TRUE)
          p = ggplot(data=d) + {plot_cmd} +
              theme_bw(base_size=10) +
              theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=10)) +
              ylab('Precision') + labs(colour = "", fill="") +
              xlab('Recall')
          ggsave("{output}", dpi=300, units="cm", height=3.5, width=10)
          """)
