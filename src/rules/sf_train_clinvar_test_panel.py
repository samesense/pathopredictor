"""Train w/ clinvar and test w/ panel"""
rule score_train_clinvar_test_panel:
    input:  DATA + 'interim/full/clinvar.dat',
            DATA + 'interim/full/panel.dat',
            DATA + 'interim/panel_genes/panel.tab'
    output: WORK + 'train_clinvar_test_panel_single{single}/roc_df_within/{cols}',
            WORK + 'train_clinvar_test_panel_single{single}/roc_df_indep/{cols}'
    shell:  'python {SCRIPTS}train_clinvar_test_panel.py {wildcards.cols} {input} {output} {wildcards.single}'

rule tmp_c:
    input: expand(WORK + 'train_clinvar_test_panel_single{single}/roc_df_within/{cols}', single=(True, False), cols=COMBO_FEATS)
