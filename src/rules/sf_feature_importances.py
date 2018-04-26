"""Compare feature importances for disease and panel vs clinvar"""

rule feature_importance:
    input:  expand(DATA + 'interim/{eval_set}.dat', eval_set=('panel', 'clinvar') )
    output: DATA + 'interim/plot_data/importances'
    shell:  'python {SCRIPTS}feature_importance.py {input} {output}'

# rule collapse_feature_importance:
#     input: expand(DATA + 'interim/importances/{eval_set}.feat_importance', eval_set=('panel',))
#     output: DATA + 'interim/plot_data/importance'
