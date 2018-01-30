"""Annotate clinvar"""
import pandas as pd

rule snpeff_clinvar:
    input:  CLINVAR
    output: DATA + 'interim/clinvar/clinvar.eff.vcf'
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {EFF} eff \
               -strict -noStats hg19 -c {EFF_CONFIG} \
               {input} > {output}"""

# rule limit_clinvar_genes:
#     input:  DATA + 'raw/EPIv6.xlsx',
#             DATA + 'interim/clinvar/clinvar.eff.vcf'
#     output: DATA + 'interim/clinvar/clinvar.use.vcf'
#     shell:  'python {SCRIPTS}limit_clinvar_genes.py {input} {output}'

rule annotateDbnsfp_clinvar:
    input:  DATA + 'interim/clinvar/clinvar.eff.vcf'
    output: DATA + 'interim/clinvar/clinvar.eff.dbnsfp.vcf'
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {SIFT} dbnsfp -v \
               -db {SIFT_DBNSFP} -f {DBNSFP_FIELDS} {input} > {output}"""

# fix pfam
# /mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/hg19.pfam.ucscgenes.enum.bed.gz
# ann fixed pfam
# parse genes
rule vcfanno_clinvar:
    input:   vcf = DATA + 'interim/clinvar/clinvar.eff.dbnsfp.vcf',
             conf = CONFIG + 'vcfanno.conf',
             lua = VCFANNO_LUA_FILE
    output:  DATA + 'interim/clinvar/clinvar.anno.vcf'
    threads: 10
    shell:   """{VCFANNO} -p {threads} -base-path {GEMINI_ANNO} -lua {input.lua} \
                {input.conf} {input.vcf} > {output}"""

rule fixHeader_clinvar:
    input:  DATA + 'interim/clinvar/clinvar.anno.vcf'
    output: DATA + 'interim/clinvar/clinvar.anno.hHack.vcf'
    shell:  'python {HEADER_HCKR} {input} {output} {HEADER_FIX}'

rule parse_clinvar_vcf:
   input:  i = DATA + 'interim/clinvar/clinvar.anno.hHack.vcf'
   output: o = DATA + 'interim/clinvar/clinvar.dat'
   run:
       with open(input.i) as f, open(output.o, 'w') as fout:
           print('chrom\tpos\tref\talt\tpfam\teff\tclinSig\taf_1kg_all\tgene\tmpc\tmtr\tProtein_Change\tconfidence\trevel\tccr', file=fout)
           for line in f:
               if line[0] != '#':
                   chrom, pos, j1, ref, alt, j2, j3, info = line.strip().split('\t')

                   ccr = '-1'
                   if 'ccr_pct' in info:
                       ccr = info.split('ccr_pct=')[1].split(';')[0]

                   mtr = '0'
                   if 'mtr=' in info:
                       mtr = info.split('mtr=')[1].split(';')[0]

                   revel = '-1'
                   if 'REVEL=' in info:
                       revel = info.split('REVEL=')[1].split(';')[0]

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

                   clin_sig = info.split('CLNSIG=')[1].split(';')[0]
                   confidence = info.split('CLNREVSTAT=')[1].split(';')[0]
                   
                   protein_change_pre = info.split('EFF=')[1].split(';')[0].split('(')[1].split('|')[3].split('/')[0]
                   protein_change = convert_protein_change(protein_change_pre)

                   eff = info.split('EFF=')[1].split(';')[0].split('(')[0]
                   gene = info.split('EFF=')[1].split(';')[0].split(',')[0].split('|')[-6]
                   ls = (chrom, pos, ref, alt, pfam, eff, clin_sig, onekg, gene, mpc, mtr, protein_change, confidence, revel, ccr)
                   print('\t'.join(ls), file=fout)

def calc_final_sig_clinvar(row):
    sig_set = set(str(row['clinSig'].split('|')))
    has_benign = '2' in sig_set or '3' in sig_set
    has_path = '4' in sig_set or '5' in sig_set
    if has_path and not has_benign:
        return 1
    if not has_path and has_benign:
        return 0
    return -1

# focus genes
rule limit_eval2_clinvar:
    input:  i = DATA + 'interim/clinvar/clinvar.dat'
    output: o = DATA + 'interim/clinvar/clinvar.limit2.dat'
    run:
        df = pd.read_csv(input.i, sep='\t')
        df.loc[:, 'clin_class'] = df.apply(calc_final_sig_clinvar, axis=1)
        crit = df.apply(lambda row: row['gene'] in FOCUS_GENES and row['clin_class'] != -1, axis=1)

        df[crit].to_csv(output.o, index=False, sep='\t')

# focus genes and missense
rule limit_eval_clinvar:
    input:  i = DATA + 'interim/clinvar/clinvar.dat'
    output: o = DATA + 'interim/clinvar/clinvar.limit.dat'
    run:
        df = pd.read_csv(input.i, sep='\t')
        df.loc[:, 'clin_class'] = df.apply(calc_final_sig_clinvar, axis=1)
        crit = df.apply(lambda row: row['gene'] in FOCUS_GENES and row['eff'] == 'missense_variant' and row['clin_class'] != -1, axis=1)
        df[crit].to_csv(output.o, index=False, sep='\t')

# missense
rule limit_eval3_clinvar:
    input:  i = DATA + 'interim/clinvar/clinvar.dat'
    output: o = DATA + 'interim/clinvar/clinvar.limit3.dat'
    run:
        df = pd.read_csv(input.i, sep='\t')
        df.loc[:, 'clin_class'] = df.apply(calc_final_sig_clinvar, axis=1)
        crit = df.apply(lambda row: row['eff'] == 'missense_variant'and row['clin_class'] != -1 and row['ccr']>-1, axis=1)
        df[crit].to_csv(output.o, index=False, sep='\t')

# singlr
# mult
# exp
# apply cutoff and higher
rule limit_clinvar_by_conf:
    input:  i=DATA + 'interim/clinvar/clinvar.limit3.dat'
    output: o=DATA + 'interim/clinvar_{conf}/clinvar_{conf}.limit3.dat'
    run:
        df = pd.read_csv(input.i, sep='\t')
        if wildcards.conf=='exp':
            crit = df.apply(lambda row: len(set(row['confidence'].split('|'))) == 1
                            and list(set(row['confidence'].split('|')))[0]==wildcards.conf, axis=1)
        elif wildcards.conf=='mult':
            crit = df.apply(lambda row: len(set(row['confidence'].split('|'))) == 1
                            and (list(set(row['confidence'].split('|')))[0]==wildcards.conf
                                 or list(set(row['confidence'].split('|')))[0]=='exp'), axis=1)
        elif wildcards.conf=='single':
            crit = df.apply(lambda row: len(set(row['confidence'].split('|'))) == 1
                            and (list(set(row['confidence'].split('|')))[0]==wildcards.conf
                                 or list(set(row['confidence'].split('|')))[0]=='exp'
                                 or list(set(row['confidence'].split('|')))[0]=='mult'), axis=1)
        else:
            i = 1/0
        df[crit].to_csv(output.o, index=False, sep='\t')

rule all_confidence:
    input: expand(DATA + 'interim/clinvar_{conf}/clinvar_{conf}.limit3.dat', conf=('exp', 'mult', 'single') )

path_color = 'f8766d'
benign_color = '00bfc4'
rule clinvar_lolly:
    input:  i = DATA + 'interim/clinvar/clinvar.limit.dat'
    output: DOCS + 'plots/clinvar/{gene}.clinvar.lolly.svg'
    run:  
        df_pre = pd.read_csv(input.i, sep='\t')
        df = df_pre[ (df_pre.gene==wildcards.gene) ]
        s = set(df[df.clin_class==0]['Protein_Change'].values)
        benign_ls = [x + '#' + benign_color for x in s if str(x) != 'nan']
        s = set(df[df.clin_class==1]['Protein_Change'].values)
        path_ls = [x + '#' + path_color for x in s if str(x) != 'nan']
        shell('~/me/bin/lollipops -domain-labels=off -o={output} -f=/home/evansj/me/fonts/arial.ttf {wildcards.gene} {path_ls} {benign_ls}')

rule all_lollies_clinvar:
    input: expand(DOCS + 'plots/clinvar/{gene}.clinvar.lolly.svg', gene=FOCUS_GENES)
