"""Snakefile for making predictions with user defined sorted and normalized vcf.
   A snakemake config file is needed. See example in configs/sm_predict_ex_config.json
   Your vcf files must be in data/interim/user_preds/
"""
include: 'const.py'
include: 'sf_ann.py'

# rule cp_user_vcf:
#     input:  DATA + 'raw/user_vcf/{name}.vcf'
#     output: DATA + 'interim/user_preds/{name}.vcf'
#     shell:  'cp {input} {output}'

rule all_predictions:
    input: expand(DATA + 'interim/user_preds/{name}.eff.dbnsfp.anno.vcf', name=config['vcf_name_ls'])
