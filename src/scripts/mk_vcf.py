import sys, pandas, twobitreader
dat_file, mutalyzer_results, twobit_file, vcf_out = sys.argv[1:]
#uc = ['Input Variant', 'Errors', 'Chromosomal Variant']
mut_df = pandas.read_csv(mutalyzer_results, sep='\t')
genome = twobitreader.TwoBitFile(twobit_file)

df_pre = pandas.read_excel(dat_file)
df = df_pre[ df_pre['Ref (single-base)'] != 'no alleles found']

def mk_var_new(row):
    if not pandas.isnull(row['Transcript']):
        if '-' == row['ref']:
            # ins
            return row['Transcript'] + ':' + row['c.'] + 'ins' + row['alt']
        elif '-' == row['alt']:
            # del
            return row['Transcript'] + ':' + row['c.'] + 'del' + row['ref']
        else:
            return row['Transcript'] + ':' + row['c.'] + row['ref'] + '>' + row['alt']
    return 'broken'

def get_pos(row):
    s = int(row['Chromosomal Variant'].split('g.')[1].split('_')[0].split('A')[0].split('G')[0].split('T')[0].split('C')[0])
    return s

def fix_ref(row):
    ref = row['Ref (single-base)']
    if ',' in ref:
        ret_ref = ref.split(',')[-1].split('>')[0]
    else:
        ret_ref = ref
    return ret_ref

def fix_alt(row):
    ref = row['Ref (single-base)']
    if ',' in ref:
        return ref.split(',')[-1].split('>')[1]
    return row['Alt (single-base)']

def mk_var(row, genome):
    if '-' == row['ref']:
        # ins
        pos = row['pos']
        append_nuc = genome['chr' + row['chrom']][pos-1:pos]
        ls  = [row['chrom'], str(row['pos']), '.',
               append_nuc, append_nuc + row['alt'] ]
    elif '-' == row['alt']:
        # del
        pos = row['pos'] - 1
        append_nuc = genome['chr' + row['chrom']][pos-1:pos]
        return [row['chrom'], str(pos), '.',
                append_nuc + row['ref'], append_nuc]
    else:
        ls  = [row['chrom'], str(row['pos']), '.',
               row['ref'], row['alt'] ]
    return ls

def mk_vcf_line(row, genome, fout):
    ls = (row['clinical_class'], row['pos_fam_count'], row['neg_fam_count'])
    info = 'CLIN_CLASS=%s;POS_FAM_COUNT=%d;NEG_FAM_COUNT=%d' % ls
    var = mk_var(row, genome)
    ls  = var + ['.', '.', info]
    print('\t'.join(ls), file=fout)

def write_vcf(df, genome, vcf_out):
    with open(vcf_out, 'w') as fout:
        print('##fileformat=VCFv4.2', file=fout)
        print('##INFO=<ID=CLIN_CLASS,Number=1,Type=String,Description="pos_fam_count">', file=fout)
        print('##INFO=<ID=POS_FAM_COUNT,Number=1,Type=Integer,Description="pos_fam_count">', file=fout)
        print('##INFO=<ID=NEG_FAM_COUNT,Number=1,Type=Integer,Description="pos_fam_count">', file=fout)
        print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', file=fout)        
        df.apply(lambda row: mk_vcf_line(row, genome, fout), axis=1)

df.loc[:, 'ref'] = df.apply(fix_ref, axis=1)
df.loc[:, 'alt'] = df.apply(fix_alt, axis=1)
df.loc[:, 'Input Variant'] = df.apply(mk_var_new, axis=1)
df.loc[:, 'chrom'] = df.apply(lambda row: row['chr'][3:], axis=1)
df.loc[:, 'clinical_class'] = df.apply(lambda row: '_'.join(row['Classification'].split()), axis=1)
m =  pandas.merge(df, mut_df, on='Input Variant', how='left').dropna()
m.loc[:, 'pos'] = m.apply(get_pos, axis=1)
uc = cols = ['chrom', 'pos', 'ref', 'alt',
             'clinical_class', 'Pos Fam Cnt', 'Neg Fam Cnt', 'Homozygous Fam Cnt',
             'Hemizygous Fam Cnt', 'Heterozygous Fam Cnt']
new_cols = {'Pos Fam Cnt':'pos_fam_count',
            'Neg Fam Cnt':'neg_fam_count',
            'Homozygous Fam Cnt':'hom_fam_count',
            'Hemizygous Fam Cnt':'hemi_fam_count',
            'Heterozygous Fam Cnt':'het_fam_count'
           }
df_final = m[uc].rename(columns=new_cols).sort_values(by=['chrom', 'pos'])
write_vcf(df_final, genome, vcf_out)
