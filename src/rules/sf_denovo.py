from const import *

rule mk_denovo_vcf:
    input:  DATA + 'raw/denovo-db.variants.v.1.5.tsv',
            DATA + 'raw/EPIv6.xlsx'
    output: DATA + 'interim/denovo/denovo.vcf'
    shell:  '  python {SCRIPTS}limit_denovo_genes.py {input} {output}'
