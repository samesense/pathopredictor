"""Test each trained classifier against the best basline result w/ fisher's test."""
import pandas as pd
import sys
#import myfisher, fisher
import scipy.stats as stats

def calc_pval(row, base):
    pval_path = stats.fisher_exact([[row['CorrectPath'], row['WrongPath']],
                                    [base['CorrectPath'], base['WrongPath']]],
                                   alternative='greater')[1]
    pval_benign = stats.fisher_exact([[row['CorrectBenign'], row['WrongBenign']],
                                      [base['CorrectBenign'], base['WrongBenign']]],
                                     alternative='greater')[1]

    return stats.combine_pvalues([pval_path, pval_benign])[1]

def calc_worst_base_pval(row, base_df):
    pvals = []
    for idx, base in base_df.iterrows():
        if base['disease'] == row['disease']:
            pvals.append( calc_pval(row, base) )
    return max(pvals)

def mk_box(row):
    if row['worst_base_pval'] < .05 and 'TRAIN' in row['st']:
        return True
    return False

def main(in_file, out_file):
    df = pd.read_csv(in_file, sep='\t')
    df_use = df#[df.st=='TRAINED_revel-ccr-is_domain']
    base_crit = df.apply(lambda row: 'BASE' in row['st']
                         and not '-' in row['st'], axis=1)
    base_df = df[base_crit]
    df_use.loc[:, 'worst_base_pval'] = df_use.apply(lambda row: calc_worst_base_pval(row, base_df), axis=1)
    df_use.loc[:, 'box'] = df_use.apply(mk_box, axis=1)
    df_use.to_csv(out_file, index=False, sep='\t')

if __name__ == "__main__":
    in_file, out_file = sys.argv[1:]
    main(in_file, out_file)
