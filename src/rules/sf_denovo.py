"""Parse variatns from denovo-db:
   http://denovo-db.gs.washington.edu/denovo-db/
"""
from const import *
from p_change import *
import pandas as pd

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

rule parse_denovo_vcf:
    input:  i = DATA + 'interim/denovo/denovo.anno.hHack.vcf'
    output: o = DATA + 'interim/denovo/denovo.dat'
    run:
        with open(input.i) as f, open(output.o, 'w') as fout:
           print('chrom\tpos\tref\talt\tpfam\teff\taf_1kg_all\tgene\tmpc\tmtr\tProtein_Change', file=fout)
           for line in f:
               if line[0] != '#':
                   chrom, pos, j1, ref, alt, j2, j3, info = line.strip().split('\t')

                   mtr = '0'
                   if 'mtr=' in info:
                       mtr = info.split('mtr=')[1].split(';')[0]

                   mpc = '0'
                   if 'mpc=' in info:
                       mpc = info.split('mpc=')[1].split(';')[0]

                   if 'pfam_domain' in info:
                       pfam = info.split('pfam_domain=')[1].split(';')[0]
                   else:
                       pfam = 'fuck'
                       
                   if 'af_1kg_all=' in info:
                       onekg = info.split('af_1kg_all=')[1].split(';')[0]
                   else:
                       onekg = '0'

                   protein_change_pre = info.split('EFF=')[1].split(';')[0].split('(')[1].split('|')[3].split('/')[0]
                   protein_change = convert_protein_change(protein_change_pre)

                   eff = info.split('EFF=')[1].split(';')[0].split('(')[0]
                   gene = info.split('EFF=')[1].split(';')[0].split(',')[0].split('|')[-6]
                   ls = (chrom, pos, ref, alt, pfam, eff, onekg, gene, mpc, mtr, protein_change)
                   print('\t'.join(ls), file=fout)

# focus genes and missense                   
rule limit_eval:                   
    input:  i = DATA + 'interim/denovo/denovo.dat'
    output: o = DATA + 'interim/denovo/denovo.limit.dat'
    run:
        df = pd.read_csv(input.i, sep='\t')
        df['y'] = 1
        crit = df.apply(lambda row: row['eff'] == 'missense_variant', axis=1)        
        df[crit].to_csv(output.o, index=False, sep='\t')

path_color = 'f8766d'
rule denovo_lolly:
    input:  i = DATA + 'interim/denovo/denovo.limit.dat'
    output: DOCS + 'plots/denovo_db/{gene}.denovo_db.lolly.svg'
    run:  
        df_pre = pd.read_csv(input.i, sep='\t')
        df = df_pre[ (df_pre.gene==wildcards.gene) ]
        s = df['Protein_Change'].values
        path_ls = [x + '#' + path_color for x in s]
        shell('~/me/bin/lollipops -domain-labels=off -o={output} -f=/home/evansj/me/fonts/arial.ttf {wildcards.gene} {path_ls}')

rule all_lollies:
    input: expand(DOCS + 'plots/denovo_db/{gene}.denovo_db.lolly.svg', gene=FOCUS_GENES)
        
