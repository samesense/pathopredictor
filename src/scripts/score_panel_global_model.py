"""Compare performance
   MPC>=2
   Train w/ all clinvar
   Train w/ clinvar for these genes
   Hold one out testing data

   Does training on clinvar/panel do better than MPC>=2?
   Does training on panel do better than training with clinvar?
"""
from collections import defaultdict
import pandas as pd
import numpy, argparse
from scipy.stats import entropy
from sklearn import linear_model, metrics, tree, svm
from sklearn.neural_network import MLPClassifier
from sklearn.externals.six import StringIO
from sklearn.preprocessing import PolynomialFeatures
from sklearn.ensemble import ExtraTreesClassifier

def eval_pred(row, col):
    if row[col] == row['y']:
        if row['y'] == 1:
            return 'CorrectPath'
        return 'CorrectBenign'
    if row['y'] == 1:
        return 'WrongPath'
    return 'WrongBenign'

def eval_mpc_raw(row):
    if row['y'] == 1:
        if row['mpc']>=float(2):
            return 'CorrectPath'
        return 'WrongPath'
    if row['mpc']>=float(2):
        return 'WrongBenign'
    return 'CorrectBenign'

def print_data_stats(disease, clinvar_df_pre, disease_df, fout):
    key_cols = ['chrom', 'pos', 'ref', 'alt']
    test_keys = {':'.join([str(x) for x in v]):True for v in disease_df[key_cols].values}

    crit = clinvar_df_pre.apply(lambda row: not ':'.join([str(row[x]) for x in key_cols]) in test_keys, axis=1)
    clinvar_df = clinvar_df_pre[crit]
    print('clinvar w/o testing data - %s: %d' % (disease, len(clinvar_df)), file=fout)

    disease_genes = set(disease_df['gene'])
    crit = clinvar_df.apply(lambda row: row['gene'] in disease_genes, axis=1)
    clinvar_df_limit_genes = clinvar_df[crit]
    print('clinvar w/o testing data for %s genes: %d' % (disease, len(clinvar_df_limit_genes)), file=fout)

    disease_panel_gene_count = len(set(disease_df['gene']))
    gg = disease_df.groupby('y').size().reset_index().rename(columns={0:'size'})
    print(disease)
    print(gg)
    benign_ex = list(gg[gg.y==0]['size'])[0]
    path_ex = list(gg[gg.y==1]['size'])[0]
    print('%s gene count: %d (%d pathogenic, %d benign)' % (disease, disease_panel_gene_count, path_ex, benign_ex), file=fout)

    return clinvar_df, clinvar_df_limit_genes

def eval_disease(disease, clinvar_df_pre, disease_df, fout_stats, fout_eval):
    clinvar_df, clinvar_df_limit_genes = print_data_stats(disease, clinvar_df_pre, disease_df, fout_stats)
    cols = ['mpc']
    
    # train clinvar
    tree_clf_clinvar = tree.DecisionTreeClassifier(max_depth=1)
    X, y = clinvar_df[cols], clinvar_df['y']
    tree_clf_clinvar.fit(X, y)

    tree_clf_clinvar_limit_genes = tree.DecisionTreeClassifier(max_depth=1)
    X, y = clinvar_df_limit_genes[cols], clinvar_df_limit_genes['y']
    tree_clf_clinvar_limit_genes.fit(X, y)

    # one gene at a time
    acc_df_ls = []
    genes = set(disease_df['gene'])

    for test_gene in genes:
        sub_train_df = disease_df[disease_df.gene != test_gene]
        tree_clf_sub = tree.DecisionTreeClassifier(max_depth=1)
        X, y = sub_train_df[cols], sub_train_df['y']
        tree_clf_sub.fit(X, y)

        test_df = disease_df[disease_df.gene == test_gene]
        X_test = test_df[cols]
        preds = tree_clf_sub.predict(X_test)
        test_df['mpc_pred'] = preds
        test_df.loc[:, 'PredictionStatusMPC'] = test_df.apply(lambda row: eval_pred(row, 'mpc_pred'), axis=1)

        preds = tree_clf_clinvar.predict(X_test)
        test_df['mpc_pred_clinvar'] = preds
        test_df.loc[:, 'PredictionStatusMPC_clinvar'] = test_df.apply(lambda row: eval_pred(row, 'mpc_pred_clinvar'), axis=1)

        preds = tree_clf_clinvar_limit_genes.predict(X_test)
        test_df['mpc_pred_clinvar_limit_genes'] = preds
        test_df.loc[:, 'PredictionStatusMPC_clinvar_limit_genes'] = test_df.apply(lambda row: eval_pred(row, 'mpc_pred_clinvar_limit_genes'), axis=1)

        # apply mpc>=2
        test_df.loc[:, 'PredictionStatusMPC>2'] = test_df.apply(eval_mpc_raw, axis=1)

        acc_df_ls.append(test_df)

    test_df = pd.concat(acc_df_ls)
    metrics_ls = ('PredictionStatusMPC', 'PredictionStatusMPC_clinvar', 'PredictionStatusMPC_clinvar_limit_genes', 'PredictionStatusMPC>2')
    for metric in metrics_ls:
        counts = test_df.groupby(metric).size().reset_index().reset_index().rename(columns={0:'size'})
        d = defaultdict(int)
        for v in list(counts.values):
            #print(v)
            _, label, count = v
            ls = (disease, 'global_' + metric, label, str(count))
            d[label] = count
            print('\t'.join(ls), file=fout_eval)
        tot_bad = d['WrongPath'] + d['WrongBenign']
        ls = (disease, 'global_' + metric, 'TotWrong', str(tot_bad))
        print('\t'.join(ls), file=fout_eval)

def main(args):
    # load clinvar
    dat_file = '../data/interim/clinvar/clinvar.limit3.dat'
    clinvar_df_pre = pd.read_csv(args.clinvar, sep='\t').rename(columns={'clin_class':'y'})

    # load genedx
    disease_genedx_df = pd.read_csv(args.gene_dx, sep='\t')
    disease_genedx_df.loc[:, 'y'] = disease_genedx_df.apply(lambda row: 1 if row['class']=='P' else 0, axis=1)
    disease_genedx_df['Disease'] = 'genedx-epi'
    
    # load other disease
    other_disease_df = pd.read_csv(args.other_disease, sep='\t')
    other_disease_df.loc[:, 'y'] = other_disease_df.apply(lambda row: 1 if row['class']=='P' else 0, axis=1)

    disease_df = pd.concat([disease_genedx_df, other_disease_df])
    diseases = set(disease_df['Disease'])

    with open(args.stats_out, 'w') as fout_stats, open(args.eval_out, 'w') as fout_eval:
        print('disease\tscore_type\teval_type\tvar_count', file=fout_eval)
        for disease in diseases:
            if str(disease) != 'nan':
                eval_disease(disease, clinvar_df_pre, disease_df[disease_df.Disease == disease],
                             fout_stats, fout_eval)
    
if __name__ == "__main__":
    desc = 'Eval global mpc cutoff'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('clinvar', 'gene_dx', 'other_disease', 'stats_out', 'eval_out')
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
