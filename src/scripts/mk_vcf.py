"""Rm duplicate genome coords.
   Favor path/benign over VUS.
   Pick largeset affected count.
"""
import pandas, argparse

new_cols = ['input_var', 'clin_class', 'pos_fam', 'neg_fam', 'hom_fam']

def process_multi(r):
    rr = r[r.clin_class=='PATHOGENIC']
    if len(rr) == 1:
        return rr[new_cols]
    
    rr = r[r.clin_class=='BENIGN']
    if len(rr) == 1:
        return rr[new_cols]
    
    r.loc[:, 'pos_sum'] = r.apply(lambda row: row['hom_fam'] + row['pos_fam'], axis=1)
    max_pos = max(r['pos_sum'].values)
    rr = r[r.pos_sum==max_pos]
    if len(rr)==1:
        return rr[new_cols]
    else:
        i = 1/0
        
def choose_one(rows):
    if len(rows) == 1:
        return rows[new_cols]
    else:
        r = rows[rows.clin_class != 'VUS']
        if not r.empty:
            if len(r) == 1:
                return r[new_cols]
            else:
                return process_multi(r)
        else:
            # only vus
            process_multi(rows)

def mk_var(row):
    ls  = [row['chrom'], str(row['pos']), '.',
           row['ref'], row['alt'] ]
    return ls

def mk_vcf_line(row, fout):
    ls = (row['input_var'], row['clin_class'],
          row['pos_fam'], row['neg_fam'], row['hom_fam'])
    info = 'INIT_VAR=%s;CLIN_CLASS=%s;POS_FAM_COUNT=%d;NEG_FAM_COUNT=%d;POS_HOM_FAM_COUNT=%d' % ls
    var = mk_var(row)
    ls  = var + ['.', '.', info]
    print('\t'.join(ls), file=fout)

def write_vcf(df, vcf_out):
    with open(vcf_out, 'w') as fout:
        print('##fileformat=VCFv4.2', file=fout)
        print('##INFO=<ID=INIT_VAR,Number=1,Type=String,Description="init_var">', file=fout)
        print('##INFO=<ID=CLIN_CLASS,Number=1,Type=String,Description="pos_fam_count">', file=fout)
        print('##INFO=<ID=POS_FAM_COUNT,Number=1,Type=Integer,Description="pos_fam_count">', file=fout)
        print('##INFO=<ID=POS_HOM_FAM_COUNT,Number=1,Type=Integer,Description="hom_fam_count">', file=fout)
        print('##INFO=<ID=NEG_FAM_COUNT,Number=1,Type=Integer,Description="pos_fam_count">', file=fout)
        print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', file=fout)
        df.apply(lambda row: mk_vcf_line(row, fout), axis=1)

def main(args):
    df = pandas.read_csv(args.dat_file, sep='\t')
    g_cols = ('chrom', 'pos', 'ref', 'alt')
    df_no_dups = df.groupby(g_cols).apply(choose_one).reset_index()
    write_vcf(df_no_dups, args.vcf_out)

if __name__ == "__main__":
    desc = 'Rm dup genome coords and convert to vcf.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('dat_file', 'vcf_out',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
