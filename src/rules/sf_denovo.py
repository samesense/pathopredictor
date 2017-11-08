from const import *

rule mk_denovo_vcf:
    input:  DATA + 'raw/denovo-db.variants.v.1.5.tsv',
            DATA + 'raw/EPIv6.xlsx'
    output: DATA + 'interim/denovo/denovo.vcf'
    shell:  '  python {SCRIPTS}limit_denovo_genes.py {input} {output}'

DBNSFP_FIELDS = 'Interpro_domain,SIFT_score,Polyphen2_HVAR_pred,RadialSVM_pred,LR_pred,Reliability_index,FATHMM_pred,MutationAssessor_pred,MutationTaster_pred,phyloP100way_vertebrate,phastCons100way_vertebrate'

rule snpeff_denovo:
    input:  DATA + 'interim/denovo/denovo.vcf' 
    output: DATA + 'interim/denovo/denovo.eff.vcf'
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {EFF} eff \
               -strict -noStats hg19 -c {EFF_CONFIG} \
               {input} > {output}"""


rule annotateDbnsfp_denovo:
    input:  DATA + 'interim/denovo/denovo.eff.vcf'
    output: DATA + 'interim/denovo/denovo.use.eff.dbnsfp.vcf'
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {SIFT} dbnsfp -v \
               -db {SIFT_DBNSFP} -f {DBNSFP_FIELDS} {input} > {output}"""

# fix pfam
# /mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/hg19.pfam.ucscgenes.enum.bed.gz
# ann fixed pfam
# parse genes
rule vcfanno_denovo:
    input:   vcf = DATA + 'interim/denovo/denovo.use.eff.dbnsfp.vcf',
             conf = CONFIG + 'vcfanno.conf',
             lua = VCFANNO_LUA_FILE
    output:  DATA + 'interim/denovo/denovo.anno.vcf'
    threads: 10
    shell:   """{VCFANNO} -p {threads} -base-path {GEMINI_ANNO} -lua {input.lua} \
                {input.conf} {input.vcf} > {output}"""

HEADER_FIX = 'eff_indel_splice,1,Flag AC,1,Integer AF,1,Float dbNSFP_FATHMM_pred,.,String dbNSFP_Interpro_domain,.,String dbNSFP_LR_pred,.,String dbNSFP_phyloP100way_vertebrate,.,String dbNSFP_phastCons100way_vertebrate,.,String dbNSFP_SIFT_score,.,String dbNSFP_Reliability_index,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_Polyphen2_HVAR_pred,.,String dbNSFP_MutationTaster_pred,.,String dbNSFP_MutationAssessor_pred,.,String'

rule fixHeader_denovo:
    input:  DATA + 'interim/denovo/denovo.anno.vcf'
    output: DATA + 'interim/denovo/denovo.anno.hHack.vcf'
    shell:  'python {HEADER_HCKR} {input} {output} {HEADER_FIX}'


