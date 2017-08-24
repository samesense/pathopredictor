"""For enriched domains, cmp path vs benign"""
import pandas, argparse

def write_dat(df, label, var_type, fout):
    path = sum(df[ (df.clin_class=='LIKLEY_PATHOGENIC') | (df.clin_class=='PATHOGENIC') ]['size'])
    ben = sum(df[ (df.clin_class=='LIKELY_BENIGN') | (df.clin_class=='BENIGN') ]['size'])
    vus = sum(df[df.clin_class=='VUS']['size'])
    path_frac, path_frac_w_vus, ben_frac_w_vus = 'NA', 'NA', 'NA'
    if path+ben+vus:
        if path+ben:
            path_frac = path/(path+ben)
        path_frac_w_vus = path/(path + ben + vus)
        ben_frac_w_vus = ben/(path + ben + vus)
    ls = [str(x) for x in (var_type, label, path, ben, vus, path_frac, path_frac_w_vus, ben_frac_w_vus)]
    print('\t'.join(ls), file=fout)

def eval_enrichment(f, var_type, fout):
    df = pandas.read_csv(f, delimiter='\t')
    sig = df.apply(lambda row: row[var_type + '_qval'] < .01 and row[var_type + '_fg_gtr'] and row['pfam'] != 'none', axis=1)
    not_sig = df.apply(lambda row: row[var_type + '_qval'] > .2 and row['pfam'] != 'none', axis=1)
    sig_exac = df.apply(lambda row: row[var_type + '_qval'] < .01 and not row[var_type + '_fg_gtr'] and row['pfam'] != 'none', axis=1)
    
    g_sig = df[sig].groupby('clin_class').size().reset_index().rename(columns={0:'size'})
    write_dat(g_sig, 'fg', var_type, fout)
    
    g_sig_exac = df[sig_exac].groupby('clin_class').size().reset_index().rename(columns={0:'size'})
    write_dat(g_sig_exac, 'exac_enriched', var_type, fout)
    
    g_ns = df[not_sig].groupby('clin_class').size().reset_index().rename(columns={0:'size'})
    write_dat(g_ns, 'not_enriched', var_type, fout)

def main(args):
    with open(args.out_file, 'w') as fout:
        fields = ('eff', 'pfam_set', 'path_count', 'benign_count', 'vus_count',
                  'path_frac_wo_vus', 'path_frac_w_vus', 'benign_frac_w_vus')
        print('\t'.join(fields), file=fout)
        eval_enrichment(args.dat_file, args.var_type, fout)

if __name__ == "__main__":
    desc = 'Pull data for report.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('dat_file', 'var_type', 'out_file',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
