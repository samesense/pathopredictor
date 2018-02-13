rule mk_sample_data:
    input:  i = WORK + 'roc_df_panel/revel-ccr'
    output: o = DATA + 'interim/webtool/sample.csv'
    run:
        keys = ['Disease', 'gene', 'chrom', 'pos', 'ref', 'alt', 'Protein_Change', 'y', 'mpc_pred'] + list(feats)
        pd.read_csv(input.i, sep='\t')[keys].rename(columns={'y':'obs_class','mpc_pred':'predicted_class'}).to_csv(output.o, index=False, sep=',')
