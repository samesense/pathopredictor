"""Annotate gnomad."""
import pandas as pd

rule snpeff_gnomad:
    input:  GNOMAD
    output: temp(DATA + 'interim/gnomad/gnomad.eff.vcf')
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {EFF} eff -dataDir {DATA}raw/snpeff/data/ \
               -strict -noStats hg19 -c {EFF_CONFIG} \
               {input} > {output}"""

# fix pfam
# /mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/hg19.pfam.ucscgenes.enum.bed.gz
# ann fixed pfam
# parse genes
rule vcfanno_gnomad:
    input:   vcf = DATA + 'interim/gnomad/gnomad.eff.vcf',
             conf = CONFIG + 'vcfanno.conf',
             lua = VCFANNO_LUA_FILE
    output:  temp(DATA + 'interim/gnomad/gnomad.anno.vcf')
    threads: 10
    shell:   """{VCFANNO} -p {threads} -base-path {GEMINI_ANNO} -lua {input.lua} \
                {input.conf} {input.vcf} > {output}"""

rule fixHeader_gnomad:
    input:  DATA + 'interim/gnomad/gnomad.anno.vcf'
    output: DATA + 'interim/gnomad/gnomad.anno.hHack.vcf'
    shell:  'python {HEADER_HCKR} {input} {output} {HEADER_FIX}'

def find_missense_eff(pos, ann, csq):
    """return eff, gene, protein_change_pre, nm"""
    for acc in (ann, csq):
        for a in acc.split(','):
            ls = a.split('|')
            if 'missense_variant' == ls[1]:
                eff, gene, protein_change, nm = ls[1], ls[3], ls[10], ls[6]
                return eff, gene, protein_change, nm
    return ls[1], ls[3], ls[10], ls[6]

def find_max_af(info, pops):
    af_ls = [-1]
    for pop in pops:
        af = info.split('AF_' + pop + '=')[1].split(';')[0]
        if af != '.':
            af_ls.append(float(af))
    return max(af_ls)

rule gnomad_rare:
    """Pull missense between .3 and 1% AF"""
    input:  i = DATA + 'interim/gnomad/gnomad.anno.hHack.vcf'
    output: o = DATA + 'interim/gnomad/gnomad.rare.dat'
    run:
        pops = ('AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS')
        with open(input.i) as f, open(output.o, 'w') as fout:
           print('chrom\tpos\tref\talt\taf\tpfam\teff\taf_1kg_all\tgene\tmpc\tmtr\tnm\tProtein_Change\trevel\tccr', file=fout)
           for line in f:
               if line[0] != '#':
                   chrom, pos, j1, ref, alt, j2, j3, info = line.strip().split('\t')
                   af = find_max_af(info, pops)
                   if True:
                       if af > 0.003:# and af < .75:
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
                           ann = info.split('ANN=')[1].split(';')[0]
                           csq = info.split('CSQ=')[1].split(';')[0]
                           eff, gene, protein_change_pre, nm = find_missense_eff(pos, ann, csq)
                           try:
                               protein_change = convert_protein_change(protein_change_pre)
                           except:
                               protein_change = 'WTF'
                           ls = (chrom, pos, ref, alt, str(af), pfam, eff, onekg, gene, mpc, mtr, nm, protein_change, revel, ccr)
                           print('\t'.join(ls), file=fout)

def limit_gnomad(input, output, low, high):
    """limit allele freqs by high and low"""
    df = pd.read_csv(input, sep='\t')
    crit = df.apply(lambda row: float(row['af']) > low and float(row['af']) < high and 'missense_variant' in row['eff'] and row['ccr']>-1, axis=1)
    df[crit].to_csv(output, index=False, sep='\t')

rule gnomad_panel:
    input:  i = DATA + 'interim/gnomad/gnomad.rare.dat'
    output: o = DATA + 'interim/gnomad/gnomad.rare.panel_{low}_{high}'
    run:
        limit_gnomad(input.i, output.o, float(wildcards.low), float(wildcards.high))

rule gnomad_clinvar:
    input:  i = DATA + 'interim/gnomad/gnomad.rare.dat'
    output: o = DATA + 'interim/gnomad/gnomad.rare.clinvar_{low}_{high}'
    run:
        limit_gnomad(input.i, output.o, float(wildcards.low), float(wildcards.high))

rule shit:
    input: DATA + 'interim/gnomad/gnomad.rare.clinvar_.05_.1'
