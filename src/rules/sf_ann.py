"""Annotate epilepsy panel vcfs"""
import csv
import pandas as pd

rule fix_mutalyzer:
    input:  DATA + 'raw/mutalyzer.{panel}'
    output: DATA + 'raw/mutalyzer.{panel}.fix'
    shell:  'cut -f 1,3 {input} > {output}'

rule fix_missing_mutalyzer:
    input:  DATA + 'raw/mutalyzer.{panel}.fix',
            DATA + 'interim/{panel}.hash'
    output: DATA + 'raw/mutalyzer.{panel}.fixMiss'
    shell:  '{PY27} {SCRIPTS}use_blat_to_fix_mutalyzer.py {input} {output}'

rule mk_dat_panel_two:
    input:  DATA + 'raw/EpilepsyVariantDataForAhmadClean_090517.xlsx',
            DATA + 'raw/mutalyzer.panel_two.fixMiss',
            '/home/evansj/me/data/ucsc/hg19.2bit'
    output: DATA + 'interim/panel_two.tab'
    shell:  'python {SCRIPTS}mk_tab_clinical_panel_two.py {input} {output}'

rule mk_dat_panel_uc:
    input:  DATA + 'interim/uc.dat',
            DATA + 'raw/mutalyzer.uc.fixMiss',
            '/home/evansj/me/data/ucsc/hg19.2bit'
    output: DATA + 'interim/uc.tab'
    shell:  'python {SCRIPTS}mk_tab_uc.py {input} {output}'

rule init_hg19_cruz:
    output: TMP + 'cruz_init_hg19'
    run:
        shell('{CRUZ_PY} {SCRIPTS}init_cruzdb.py hg19')
        shell('touch {output}')

def fix_transcript(row):
    return row['Transcript'].split('.')[0]

# uc data is missing chrom
rule add_uc_chroms:
    input:  x = DATA + 'raw/UC_all_panel_variants_01_20_2016.xlsx',
            d = TMP + 'cruz_init_hg19'
    output: o = DATA + 'interim/uc.dat'
    run:
        df = pd.read_excel(input.x)
        df.loc[:, 'simple_nm_tmp'] = df.apply(fix_transcript, axis=1)
        with open(input.x + '.nm', 'w') as fout:
            for nm in set(df['simple_nm_tmp']):
                print(nm, file=fout)
        shell('{CRUZ_PY} {SCRIPTS}add_chr_to_mutalyzer.py {input.x}.nm {output}.tmp')
        df_chrom = pd.read_csv(output.o + '.tmp', sep='\t')
        m = pd.merge(df, df_chrom, how='left', on='simple_nm_tmp')
        shell('rm {input.x}.nm {output}.tmp')
        m.to_csv(output.o, sep='\t', index=False)

# rule mk_dat_panel_one:
#     input:  DATA + 'raw/EPIv6.xlsx',
#             DATA + 'raw/mut.fix',
#             '/home/evansj/me/data/ucsc/hg19.2bit'
#     output: DATA + 'interim/EPIv6.tab'
#     shell:  'python {SCRIPTS}mk_tab.py {input} {output}'

rule mk_vcf_epi:
    input:  DATA + 'interim/{lab}.tab'
    output: DATA + 'interim/{lab}.pre.vcf'
    shell:  'python {SCRIPTS}mk_vcf.py {input} {output}'

rule sort_vcf_epi:
    input:  DATA + 'interim/{lab}.pre.vcf'
    output: DATA + 'interim/epi/{lab}.vcf'
    shell:  'cat {input} | vcf-sort > {output}'

rule bgzipVcf_epi:
    input:  DATA + 'interim/epi/{lab}.vcf'
    output: DATA + 'interim/epi/{lab}.vcf.gz'
    shell:  '{BGZ} -c {input} > {output}'

rule snpeff_epi:
    input:  DATA + 'interim/epi/{lab}.vcf.gz'
    output: DATA + 'interim/epi/{lab}.eff.vcf'
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {EFF} eff \
               -strict -noStats hg19 -c {EFF_CONFIG} \
               {input} > {output}"""

rule annotateDbnsfp_epi:
    input:  DATA + 'interim/epi/{lab}.eff.vcf'
    output: DATA + 'interim/epi/{lab}.eff.dbnsfp.vcf'
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {SIFT} dbnsfp -v \
               -db {SIFT_DBNSFP} -f {DBNSFP_FIELDS} {input} > {output}"""

# fix pfam
# /mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/hg19.pfam.ucscgenes.enum.bed.gz
# ann fixed pfam
# parse genes
rule vcfanno_epi:
    input:   vcf = DATA + 'interim/epi/{lab}.eff.dbnsfp.vcf',
             conf = CONFIG + 'vcfanno.conf',
             lua = VCFANNO_LUA_FILE
    output:  DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.vcf'
    threads: 10
    shell:   """{VCFANNO} -p {threads} -base-path {GEMINI_ANNO} -lua {input.lua} \
                {input.conf} {input.vcf} > {output}"""

rule fixHeader_epi:
    input:  DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.vcf'
    output: DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.hHack.vcf'
    shell:  'python {HEADER_HCKR} {input} {output} {HEADER_FIX}'

rule zip_epi:
    input:  DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.hHack.vcf'
    output: DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.hHack.vcf.gz'
    shell:  '{BGZ} -c {input} > {output}'

# neg fam counts ppl
# pos fam counts ppl
# need to convert hom to a count of two
rule parse_vcf_epi:
   input:  i = DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.hHack.vcf'
   output: o = DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.hHack.dat.xls',
           o2 = DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.hHack.splitPfam.dat'
   run:
       with open(input.i) as f, open(output.o, 'w') as fout, open(output.o2, 'w') as fout_split_pfam:
           print('chrom\tpos\tref\talt\tclin_class\tpfam\taf_1kg_all\teff\tpos_fam\tneg_fam\tgene\tmpc\tmtr\trevel\texac_af\texac_ac\texac_an\texac_cov_frac\tkaviar_af\tc.\tProtein_Change\tHugo_Symbol\tccr', file=fout)
           print('chrom\tpos\tref\talt\tclin_class\tpfam\taf_1kg_all\teff\tpos_fam\tneg_fam\tgene\tmpc\tmtr\trevel\tccr',
                 file=fout_split_pfam)
           for line in f:
               if not line[0] == '#':
                   chrom, pos, j1, ref, alt, j2, j3, info = line.strip().split('\t')

                   c_dot = info.split('INIT_VAR=')[1].split(';')[0]

                   exac_af = '0'
                   kv_af = '0'
                   exac_cov_frac = '0'
                   exac_ac = '0'
                   exac_an = '0'
                   if 'af_exac_all=' in info:
                       exac_af = info.split('af_exac_all=')[1].split(';')[0]
                   if 'kv_af=' in info:
                       kv_af = info.split('kv_af=')[1].split(';')[0]
                   if 'totExacCov_10=' in info:
                       exac_cov_frac = info.split('totExacCov_10=')[1].split(';')[0]
                   if 'an_exac_all' in info:
                       exac_an = info.split('an_exac_all=')[1].split(';')[0]
                   if 'ac_exac_all' in info:
                       exac_ac = info.split('ac_exac_all=')[1].split(';')[0]

                   ccr = '-1'
                   if 'ccr_pct' in info:
                       ccr = info.split('ccr_pct=')[1].split(';')[0]

                   mpc = '0'
                   if 'mpc=' in info:
                       mpc = info.split('mpc=')[1].split(';')[0]

                   mtr = '0'
                   if 'mtr=' in info:
                       mtr = info.split('mtr=')[1].split(';')[0]

                   revel = '-1'
                   if 'REVEL=' in info:
                       revel = info.split('REVEL=')[1].split(';')[0]


                   clin = info.split('CLIN_CLASS=')[1].split(';')[0]

                   pos_fam = int(info.split('POS_FAM_COUNT=')[1].split(';')[0])
                   neg_fam = info.split('NEG_FAM_COUNT=')[1].split(';')[0]
                   hom_fam = int(info.split('POS_HOM_FAM_COUNT=')[1].split(';')[0])
                   pos_fam = str(pos_fam + hom_fam)

                   if 'pfam_domain' in info:
                       pfam = info.split('pfam_domain=')[1].split(';')[0]
                   else:
                       pfam = 'fuck'

                   if 'af_1kg_all=' in info:
                       onekg = info.split('af_1kg_all=')[1].split(';')[0]
                   else:
                       onekg = '0'

                   eff = info.split('EFF=')[1].split(';')[0].split('(')[0]
                   # (MODERATE|MISSENSE|Ata/Gta|p.Ile199Val/
                   #protein_change = 'NA'
                   #if eff == 'missense_variant':
                   protein_change_pre = info.split('EFF=')[1].split(';')[0].split('(')[1].split('|')[3].split('/')[0]
                   protein_change = convert_protein_change(protein_change_pre)

                   gene = info.split('EFF=')[1].split(';')[0].split(',')[0].split('|')[-6]
                   ls = (chrom, pos, ref, alt, clin, pfam, onekg, eff, pos_fam,
                         neg_fam, gene, mpc, mtr, revel, exac_af, exac_ac, exac_an,
                         exac_cov_frac, kv_af, c_dot, protein_change, gene, ccr)
                   print('\t'.join(ls), file=fout)

                   for p in pfam.split(','):
                       ls = (chrom, pos, ref, alt, clin, p, onekg, eff, pos_fam, neg_fam, gene, mpc, mtr, revel, ccr)
                       print('\t'.join(ls), file=fout_split_pfam)

def mk_class_epi(row):
    if row['clin_class'] in ('Benign', 'BENIGN', 'LIKELY BENIGN', 'likely benign', 'benign', 'LIKELY_BENIGN', 'likely_benign'):
        return 'B'
    elif row['clin_class'] in ('pathogenic recessive', 'pathogenic dominant', 'likely pathogenic', 'LIKLEY_PATHOGENIC', 'pathogenic_recessive',
                               'LIKELY_PATHOGENIC', 'likely_pathogenic', 'likely_pathogenic', 'Reduced_function_allele', 'pathogenic_dominant',
                               'PATHOGENIC', 'LIKELY PATHOGENIC', 'Reduced function allele', 'pathogenic'):
        return 'P'
    elif str(row['clin_class']) in ('VUS', 'nan', 'VOUS'):
        return 'V'
    else:
        print(str(row['clin_class']))
        i = 1/0

# focus genes and missense
# must have ccr score
rule limit_eval_epi:
    input:  i = DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.hHack.dat.xls'
    output: o = DATA + 'interim/epi/{lab}.eff.dbnsfp.anno.hHack.dat.limit.xls'
    run:
        df = pd.read_csv(input.i, sep='\t')
        df.loc[:, 'class'] = df.apply(mk_class_epi, axis=1)
        crit = df.apply(lambda row: row['eff'] == 'missense_variant' and row['class'] != 'V' and row['ccr']>-1, axis=1)
        df[crit].to_csv(output.o, index=False, sep='\t')

path_color = 'f8766d'
benign_color = '00bfc4'
# rule denovo_lolly_epi:
#     input:  i = DATA + 'interim/{lab}.eff.dbnsfp.anno.hHack.dat.limit.xls'
#     output: DOCS + 'plots/{lab}/{gene}.{lab}.lolly.svg'
#     run:
#         df_pre = pd.read_csv(input.i, sep='\t')
#         df = df_pre[ (df_pre.gene==wildcards.gene) ]
#         s = set(df[df['class']=='B']['Protein_Change'].values)
#         benign_ls = [x + '#' + benign_color for x in s]
#         s = set(df[df['class']=='P']['Protein_Change'].values)
#         path_ls = [x + '#' + path_color for x in s]
#         shell('~/me/bin/lollipops -domain-labels=off -o={output} -f=/home/evansj/me/fonts/arial.ttf {wildcards.gene} {path_ls} {benign_ls}')

# rule all_lollies_epi:
#     input: expand(DOCS + 'plots/{panel}/{gene}.{panel}.lolly.svg', gene=FOCUS_GENES, panel=('EPIv6', 'panel_two', 'uc'))

# rule all_labs_epi:
#     input: expand(DATA + 'interim/{lab}.eff.dbnsfp.anno.hHack.dat.limit.xls', lab=('EPIv6', 'uc',))
