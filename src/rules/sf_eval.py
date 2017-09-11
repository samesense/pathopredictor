"""Make predictions"""

rule clinvar_eval:
    input:  DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.xls',
            DATA + 'interim/clinvar/clinvar.dat'
    output: DOCS + 'plots/missense_clinvar_roc_feature_union.png',
            WORK + 'eval/auc/missense_clinvar_roc_feature_union',
            WORK + 'eval/missense_clinvar_roc_feature_union.dat'
    shell:  'python {SCRIPTS}clinvar_eval.py {input} {output}'
