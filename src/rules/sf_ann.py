"""Annotate vcf"""
from const import *

rule mk_dat:
    input:  DATA + 'raw/EPIv6.xlsx',
            DATA + 'raw/mut.fix',
            '/home/evansj/me/data/ucsc/hg19.2bit'
    output: DATA + 'interim/EPIv6.tab'
    shell:  'python {SCRIPTS}mk_tab.py {input} {output}'

rule mk_vcf:
    input:  DATA + 'interim/EPIv6.tab'
    output: DATA + 'interim/EPIv6.pre.vcf'
    shell:  'python {SCRIPTS}mk_vcf.py {input} {output}'

rule sort_vcf:
    input:  DATA + 'interim/EPIv6.pre.vcf'
    output: DATA + 'interim/EPIv6.vcf'
    shell:  'cat {input} | vcf-sort > {output}'

rule bgzipVcf:
    input:  DATA + 'interim/EPIv6.vcf'
    output: DATA + 'interim/EPIv6.vcf.gz'
    shell:  '{BGZ} -c {input} > {output}'

rule snpeff:
    input:  DATA + 'interim/EPIv6.vcf.gz'
    output: DATA + 'interim/EPIv6.eff.vcf'
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {EFF} eff \
               -strict -noStats hg19 -c {EFF_CONFIG} \
               {input} > {output}"""

DBNSFP_FIELDS = 'Interpro_domain,SIFT_score,Polyphen2_HVAR_pred,RadialSVM_pred,LR_pred,Reliability_index,FATHMM_pred,MutationAssessor_pred,MutationTaster_pred,phyloP100way_vertebrate,phastCons100way_vertebrate'

rule annotateDbnsfp:
    input:  DATA + 'interim/EPIv6.eff.vcf'
    output: DATA + 'interim/EPIv6.eff.dbnsfp.vcf'
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {SIFT} dbnsfp -v \
               -db {SIFT_DBNSFP} -f {DBNSFP_FIELDS} {input} > {output}"""

# fix pfam
# /mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/hg19.pfam.ucscgenes.enum.bed.gz
# ann fixed pfam
# parse genes
rule vcfanno:
    input:   vcf = DATA + 'interim/EPIv6.eff.dbnsfp.vcf',
             conf = CONFIG + 'vcfanno.conf',
             lua = VCFANNO_LUA_FILE
    output:  DATA + 'interim/EPIv6.eff.dbnsfp.anno.vcf'
    threads: 10
    shell:   """{VCFANNO} -p {threads} -base-path {GEMINI_ANNO} -lua {input.lua} \
                {input.conf} {input.vcf} > {output}"""

HEADER_FIX = 'eff_indel_splice,1,Flag AC,1,Integer AF,1,Float dbNSFP_FATHMM_pred,.,String dbNSFP_Interpro_domain,.,String dbNSFP_LR_pred,.,String dbNSFP_phyloP100way_vertebrate,.,String dbNSFP_phastCons100way_vertebrate,.,String dbNSFP_SIFT_score,.,String dbNSFP_Reliability_index,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_Polyphen2_HVAR_pred,.,String dbNSFP_MutationTaster_pred,.,String dbNSFP_MutationAssessor_pred,.,String'

rule fixHeader:
    input:  DATA + 'interim/EPIv6.eff.dbnsfp.anno.vcf'
    output: DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.vcf'
    shell:  'python {HEADER_HCKR} {input} {output} {HEADER_FIX}'

rule zip:
    input:  DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.vcf'
    output: DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.vcf.gz'
    shell:  '{BGZ} -c {input} > {output}'

# rule gemini:
#     input:  df=DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.vcf.gz'
#             jp=PWD + 'files/CAP.JUNK_PED.ped'
#     output: DATA + 'interim/EPIv6.db'
#     shell:  """{GPY} {VCFTODB} --legacy-compression {input.df} {input.jp} {output}"""

# neg fam counts ppl
# pos fam counts ppl
# need to convert hom to a count of two
rule parse_vcf:
   input:  i = DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.vcf'
   output: o = DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat',
           o2 = DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.dat'
   run:
       with open(input.i) as f, open(output.o, 'w') as fout, open(output.o2, 'w') as fout_split_pfam:
           print('chrom\tpos\tref\talt\tclin_class\tpfam\taf_1kg_all\teff\tgene\tmpc', file=fout)
           print('chrom\tpos\tref\talt\tclin_class\tpfam\taf_1kg_all\teff\tpos_fam\tneg_fam\tgene\tmpc',
                 file=fout_split_pfam)
           for line in f:
               if not line[0] == '#':
                   chrom, pos, j1, ref, alt, j2, j3, info = line.strip().split('\t')

                   mpc = 'NA'
                   if 'mpc=' in info:
                       mpc = info.split('mpc=')[1].split(';')[0]

                   clin = info.split('CLIN_CLASS=')[1].split(';')[0]

                   pos_fam = int(info.split('POS_FAM_COUNT=')[1].split(';')[0])
                   neg_fam = info.split('NEG_FAM_COUNT=')[1].split(';')[0]
                   hom_fam = int(info.split('POS_HOM_FAM_COUNT=')[1].split(';')[0])
                   pos_fam = str(pos_fam + hom_fam)

                   if 'pfam_domain' in info:
                       pfam = info.split('pfam_domain=')[1].split(';')[0]
                   else:
                       pfam = 'none'

                   if 'af_1kg_all=' in info:
                       onekg = info.split('af_1kg_all=')[1].split(';')[0]
                   else:
                       onekg = '0'

                   eff = info.split('EFF=')[1].split(';')[0].split('(')[0]
                   gene = info.split('EFF=')[1].split(';')[0].split(',')[0].split('|')[-6]
                   ls = (chrom, pos, ref, alt, clin, pfam, onekg, eff, pos_fam, neg_fam, gene, mpc)
                   print('\t'.join(ls), file=fout)

                   for p in pfam.split(','):
                       ls = (chrom, pos, ref, alt, clin, p, onekg, eff, pos_fam, neg_fam, gene, mpc)
                       print('\t'.join(ls), file=fout_split_pfam)
                   
