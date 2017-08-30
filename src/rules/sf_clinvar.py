"""Annotate clinvar"""
from const import *

rule snpeff_clinvar:
    input:  CLINVAR
    output: DATA + 'interim/clinvar/clinvar.eff.vcf'
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {EFF} eff \
               -strict -noStats hg19 -c {EFF_CONFIG} \
               {input} > {output}"""

rule limit_clinvar_genes:
    input:  DATA + 'raw/EPIv6.xlsx',
            DATA + 'interim/clinvar/clinvar.eff.vcf'
    output: DATA + 'interim/clinvar/clinvar.use.vcf'
    shell:  'python {SCRIPTS}limit_clinvar_genes.py {input} {output}'

DBNSFP_FIELDS = 'Interpro_domain,SIFT_score,Polyphen2_HVAR_pred,RadialSVM_pred,LR_pred,Reliability_index,FATHMM_pred,MutationAssessor_pred,MutationTaster_pred,phyloP100way_vertebrate,phastCons100way_vertebrate'

rule annotateDbnsfp_clinvar:
    input:  DATA + 'interim/clinvar/clinvar.use.vcf'
    output: DATA + 'interim/clivnar/clinvar.use.eff.dbnsfp.vcf'
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {SIFT} dbnsfp -v \
               -db {SIFT_DBNSFP} -f {DBNSFP_FIELDS} {input} > {output}"""

# fix pfam
# /mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/hg19.pfam.ucscgenes.enum.bed.gz
# ann fixed pfam
# parse genes
rule vcfanno_clinvar:
    input:   vcf = DATA + 'interim/clivnar/clinvar.use.eff.dbnsfp.vcf',
             conf = CONFIG + 'vcfanno.conf',
             lua = VCFANNO_LUA_FILE
    output:  DATA + 'interim/clinvar/clinvar.anno.vcf'
    threads: 10
    shell:   """{VCFANNO} -p {threads} -base-path {GEMINI_ANNO} -lua {input.lua} \
                {input.conf} {input.vcf} > {output}"""

HEADER_FIX = 'eff_indel_splice,1,Flag AC,1,Integer AF,1,Float dbNSFP_FATHMM_pred,.,String dbNSFP_Interpro_domain,.,String dbNSFP_LR_pred,.,String dbNSFP_phyloP100way_vertebrate,.,String dbNSFP_phastCons100way_vertebrate,.,String dbNSFP_SIFT_score,.,String dbNSFP_Reliability_index,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_Polyphen2_HVAR_pred,.,String dbNSFP_MutationTaster_pred,.,String dbNSFP_MutationAssessor_pred,.,String'

rule fixHeader_clinvar:
    input:  DATA + 'interim/clinvar/clinvar.anno.vcf'
    output: DATA + 'interim/clinvar/clinvar.anno.hHack.vcf'
    shell:  'python {HEADER_HCKR} {input} {output} {HEADER_FIX}'

rule parse_clinvar_vcf:
   input:  i = DATA + 'interim/clinvar/clinvar.anno.hHack.vcf'
   output: o = DATA + 'interim/clinvar/clinvar.dat'
   run:
       with open(input.i) as f, open(output.o, 'w') as fout:
           print('chrom\tpos\tref\talt\tpfam\teff\tclinSig\taf_1kg_all\tgene\tmpc', file=fout)
           for line in f:
               if line[0] != '#':
                   chrom, pos, j1, ref, alt, j2, j3, info = line.strip().split('\t')

                   mpc = 'NA'
                   if 'mpc=' in info:
                       mpc = info.split('mpc=')[1].split(';')[0]

                   if 'pfam_domain' in info:
                       pfam = info.split('pfam_domain=')[1].split(';')[0]
                   else:
                       pfam = 'none'
                       
                   if 'af_1kg_all=' in info:
                       onekg = info.split('af_1kg_all=')[1].split(';')[0]
                   else:
                       onekg = '0'

#                   print(info)
                   clin_sig = info.split('CLNSIG=')[1].split(';')[0]

                   eff = info.split('EFF=')[1].split(';')[0].split('(')[0]
                   gene = info.split('EFF=')[1].split(';')[0].split(',')[0].split('|')[-6]
                   ls = (chrom, pos, ref, alt, pfam, eff, clin_sig, onekg, gene, mpc)
                   print('\t'.join(ls), file=fout)

