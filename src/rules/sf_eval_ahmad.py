"""Combine evaluations by disease."""
import pandas as pd

rule percent_wrong:
    input:  panel = WORK + '{method}.eval_panel.eval',
            non_panel = expand( WORK + '{{method}}.eval_{dat}', dat=('clinvar', 'denovo', 'clinvar_mult', 'clinvar_single', 'clinvar_exp') )
    output: o = WORK + '{method}.eval_panel.eval.percentWrong'
    run:
        def calc_wrong_percent(rows):
            tot_wrong = rows[rows.eval_type == 'TotWrong']['var_count']
            tot_preds = sum(rows[rows.eval_type != 'TotWrong']['var_count'])
            return tot_wrong/tot_preds

        def load_other_df(afile):
            panel_df = (pd.read_csv(input.panel, sep='\t')
                        .groupby(('clinvar_type', 'disease', 'score_type'))
                        .apply(calc_wrong_percent)
                        .reset_index()
                        .rename(columns={0:'percent_wrong'})
                        .to_csv(output.o, index=False, sep='\t') )
            return panel_df

        non_panel_dfs = pd.concat( [load_other_df(afile) for afile in input.non_panel] )
        panel_df = (pd.read_csv(input.panel, sep='\t')
                    .groupby(('disease', 'score_type'))
                    .apply(calc_wrong_percent)
                    .reset_index()
                    .rename(columns={0:'percent_wrong'})
                    .to_csv(output.o, index=False, sep='\t') )
        pd.merge(panel_df, non_panel_dfs, on='disease', how='left').to_csv(output.o, index=False, sep='\t')

                   
    
