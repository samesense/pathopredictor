"""For enriched domains, cmp path vs benign"""
import pandas, argparse

def eval_enrichment(f, var_type, fout):
    df = pandas.read_csv(f, delimiter='\t')
    sig = df.apply(lambda row: row[var_type + '_qval'] < .01 and row[var_type + '_fg_gtr'] and row['pfam'] != 'none', axis=1)
    not_sig = df.apply(lambda row: row[var_type + '_qval'] > .2 and row['pfam'] != 'none', axis=1)
    sig_exac = df.apply(lambda row: row[var_type + '_qval'] < .01 and not row[var_type + '_fg_gtr'] and row['pfam'] != 'none', axis=1)
    
    g_sig = df[sig].groupby('clin_class').size().reset_index().rename(columns={0:'size'})
    path = sum(g_sig[ (g_sig.clin_class=='LIKELY_PATHOGENIC') | (g_sig.clin_class=='PATHOGENIC') ]['size'])
    ben = sum(g_sig[ (g_sig.clin_class=='LIKELY_BENIGN') | (g_sig.clin_class=='BENIGN') ]['size'])
    frac = 'NA'
    if path+ben:
        frac = path/(path+ben)
    ls = [str(x) for x in (var_type, 'fg', path, ben, frac)]
    print('\t'.join(ls), file=fout)
    
    g_sig_exac = df[sig_exac].groupby('clin_class').size().reset_index().rename(columns={0:'size'})
    path = sum(g_sig_exac[ (g_sig_exac.clin_class=='LIKELY_PATHOGENIC') | (g_sig_exac.clin_class=='PATHOGENIC') ]['size'])
    ben = sum(g_sig_exac[ (g_sig_exac.clin_class=='LIKELY_BENIGN') | (g_sig_exac.clin_class=='BENIGN') ]['size'])
    frac = 'NA'
    if path+ben:
        frac = path/(path+ben)
    ls = [str(x) for x in (var_type, 'exac_enriched', path, ben, frac)]
    print('\t'.join(ls), file=fout)
    
    g_ns = df[not_sig].groupby('clin_class').size().reset_index().rename(columns={0:'size'})
    path = sum(g_ns[ (g_ns.clin_class=='LIKELY_PATHOGENIC') | (g_ns.clin_class=='PATHOGENIC') ]['size'])
    ben = sum(g_ns[ (g_ns.clin_class=='LIKELY_BENIGN') | (g_ns.clin_class=='BENIGN') ]['size'])
    frac = 'NA'
    if path+ben:
        frac = path/(path+ben)
    ls = [str(x) for x in (var_type, 'not_enriched', path, ben, frac)]
    print('\t'.join(ls), file=fout)

def main(args):
    with open(args.out_file, 'w') as fout:
        print('eff\tpfam_set\tpath_count\tbenign_count\tpath_frac', file=fout)
        eval_enrichment(args.dat_file, args.var_type, fout)

if __name__ == "__main__":
    desc = 'Pull data for report.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('dat_file', 'var_type', 'out_file',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
