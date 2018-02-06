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

def mk_basic_call(row, score_cols):
    """Call is pathogenic/1 if all scores are met"""
    cutoffs = {'is_domain':1, 'mpc':2, 'revel':.375, 'ccr':90}
    for col in score_cols:
        if row[col] < cutoffs[col]:
            return 0
    return 1

def eval_mpc_raw(row, score_cols):
    call = mk_basic_call(row, score_cols)
    if row['y'] == 1:
        if call == 1:
            return 'CorrectPath'
        return 'WrongPath'
    # else benign
    if call == 1:
        return 'WrongBenign'
    return 'CorrectBenign'

def clinvar_stats(disease, clinvar_df_pre, disease_df, fout, clinvar_label):
    key_cols = ['chrom', 'pos', 'ref', 'alt']
    test_keys = {':'.join([str(x) for x in v]):True for v in disease_df[key_cols].values}

    crit = clinvar_df_pre.apply(lambda row: not ':'.join([str(row[x]) for x in key_cols]) in test_keys, axis=1)
    clinvar_df = clinvar_df_pre[crit]
    print('%s clinvar w/o testing data - %s: %d' % (clinvar_label, disease, len(clinvar_df)), file=fout)

    disease_genes = set(disease_df['gene'])
    crit = clinvar_df.apply(lambda row: row['gene'] in disease_genes, axis=1)
    clinvar_df_limit_genes = clinvar_df[crit]
    print('%s clinvar w/o testing data for %s genes: %d' % (clinvar_label, disease, len(clinvar_df_limit_genes)), file=fout)

    return clinvar_df, clinvar_df_limit_genes

def print_data_stats(disease, clinvar_df_pre_ls, disease_df, fout, clin_labels):
    clin_dat = list(map(lambda x: clinvar_stats(disease, x[0], disease_df, fout, x[1]),
                        list(zip(clinvar_df_pre_ls, clin_labels))))

    disease_panel_gene_count = len(set(disease_df['gene']))
    gg = disease_df.groupby('y').size().reset_index().rename(columns={0:'size'})
    benign_ex = list(gg[gg.y==0]['size'])[0]
    path_ex = list(gg[gg.y==1]['size'])[0]
    print('%s gene count: %d (%d pathogenic, %d benign)' % (disease, disease_panel_gene_count, path_ex, benign_ex), file=fout)

    return [x[0] for x in clin_dat], [x[1] for x in clin_dat]

def eval_clinvar(label, cols, clinvar_df, disease_df):
    #print(label, clinvar_df.columns.values)
    if len(clinvar_df):
        # train clinvar
        tree_clf_clinvar = tree.DecisionTreeClassifier( max_depth=len(cols) )
        X, y = clinvar_df[cols], clinvar_df['y']
        tree_clf_clinvar.fit(X, y)

        X_test = disease_df[cols]
        print(X_test)
        preds = tree_clf_clinvar.predict(X_test)

        disease_df['mpc_pred_clinvar_' + label] = preds
        disease_df.loc[:, 'PredictionStatusMPC_clinvar_' + label] = disease_df.apply(lambda row: eval_pred(row, 'mpc_pred_clinvar_' + label), axis=1)

def print_eval(disease, test_df, metric, fout):
    if metric in test_df.columns.values:
        counts = test_df.groupby(metric).size().reset_index().reset_index().rename(columns={0:'size'})
        d = defaultdict(int)
        for v in list(counts.values):
            _, label, count = v
            ls = (disease, 'global_' + metric, label, str(count))
            d[label] = count
            print('\t'.join(ls), file=fout)
        tot_bad = d['WrongPath'] + d['WrongBenign']
        ls = (disease, 'global_' + metric, 'TotWrong', str(tot_bad))
        print('\t'.join(ls), file=fout)

def eval_disease(disease, clinvar_df_pre_ls, disease_df, fout_stats, fout_eval, clin_labels, cols):
    #plot_panel_tree(disease, disease_df, cols)
    clinvar_df_ls, clinvar_df_limit_genes_ls = print_data_stats(disease, clinvar_df_pre_ls, disease_df, fout_stats, clin_labels)
    list(map(lambda x: eval_clinvar(x[0], cols, x[1], disease_df),
             list(zip(clin_labels, clinvar_df_ls))))
    list(map(lambda x: eval_clinvar(x[0], cols, x[1], disease_df),
             list(zip([x + '_limitGene' for x in clin_labels], clinvar_df_limit_genes_ls))))

    # apply mpc>=2
    disease_df.loc[:, 'PredictionStatusMPC>2'] = disease_df.apply(lambda x: eval_mpc_raw(x, cols), axis=1)

    # one gene at a time
    acc_df_ls = []
    genes = set(disease_df['gene'])

    metric_ls = ['PredictionStatusMPC_clinvar_' + label for label in clin_labels] + ['PredictionStatusMPC_clinvar_' + label + '_limitGene' for label in clin_labels] + ['PredictionStatusMPC>2']
    list(map(lambda x: print_eval(disease, disease_df, x, fout_eval), metric_ls))

    for test_gene in genes:
        sub_train_df = disease_df[disease_df.gene != test_gene]
        tree_clf_sub = tree.DecisionTreeClassifier( max_depth=len(cols) )
        X, y = sub_train_df[cols], sub_train_df['y']
        tree_clf_sub.fit(X, y)

        #regression for multiple scores
        if len(cols)>1:
            lm =  linear_model.LinearRegression(normalize=True, fit_intercept=True)
            lm.fit(X, y)

        test_df = disease_df[disease_df.gene == test_gene]
        X_test = test_df[cols]
        preds = tree_clf_sub.predict(X_test)
        if len(cols)>1:
            lm_preds = lm.predict(X_test)
            test_df['-'.join(cols) + '_pred_lm'] = lm_preds
        test_df['mpc_pred'] = preds
        test_df.loc[:, 'PredictionStatusMPC'] = test_df.apply(lambda row: eval_pred(row, 'mpc_pred'), axis=1)

        acc_df_ls.append(test_df)

    test_df = pd.concat(acc_df_ls)
    print_eval(disease, test_df, 'PredictionStatusMPC', fout_eval)
    return test_df

def main(args):
    score_cols = args.score_cols.split('-')
    FOCUS_GENES = ('SCN1A','SCN2A','KCNQ2', 'KCNQ3', 'CDKL5',
                   'PCDH19', 'SCN1B', 'SCN8A', 'SLC2A1',
                   'SPTAN1', 'STXBP1', 'TSC1')

    # load clinvar
    clin_labels = ('tot', 'single', 'mult', 'exp', 'denovo')
    clinvars = (args.clinvar, args.clinvar_single, args.clinvar_mult, args.clinvar_exp)
    clinvar_df_pre_ls = list(map(lambda x: pd.read_csv(x, sep='\t').rename(columns={'clin_class':'y'}), clinvars))
    # add denvo
    clinvar_df_pre_ls.append( pd.read_csv(args.denovo, sep='\t') )

    for df in clinvar_df_pre_ls:
        df.loc[:, 'is_domain'] = df.apply(lambda row: 0 if 'none' in row else 1, axis=1)
    # load genedx epi
    disease_genedx_df = pd.read_csv(args.gene_dx, sep='\t')
    disease_genedx_df.loc[:, 'y'] = disease_genedx_df.apply(lambda row: 1 if row['class']=='P' else 0, axis=1)
    disease_genedx_df['Disease'] = 'genedx-epi'
    crit = disease_genedx_df.apply(lambda row: row['gene'] in FOCUS_GENES, axis=1)
    disease_genedx_limitGene_df= disease_genedx_df[crit]
    disease_genedx_limitGene_df['Disease'] = 'genedx-epi-limitGene'

    # load uc epi
    disease_uc_df = pd.read_csv(args.uc, sep='\t')
    disease_uc_df.loc[:, 'y'] = disease_uc_df.apply(lambda row: 1 if row['class']=='P' else 0, axis=1)
    disease_uc_df['Disease'] = 'uc-epi'
    crit = disease_uc_df.apply(lambda row: row['gene'] in FOCUS_GENES, axis=1)
    disease_uc_limitGene_df= disease_uc_df[crit]
    disease_uc_limitGene_df['Disease'] = 'uc-epi-limitGene'

    # load other disease
    other_disease_df = pd.read_csv(args.other_disease, sep='\t')
    other_disease_df.loc[:, 'y'] = other_disease_df.apply(lambda row: 1 if row['class']=='P' else 0, axis=1)

    disease_df = pd.concat([disease_genedx_df, disease_uc_df,
                            disease_genedx_limitGene_df, disease_uc_limitGene_df,
                            other_disease_df])
    disease_df.loc[:, 'is_domain'] = disease_df.apply(lambda row: 0 if 'none' in row['pfam'] else 1, axis=1)
    diseases = set(disease_df['Disease'])

    eval_df_ls = []
    with open(args.stats_out, 'w') as fout_stats, open(args.eval_out, 'w') as fout_eval:
        print('disease\tscore_type\teval_type\tvar_count', file=fout_eval)
        for disease in diseases:
            if str(disease) != 'nan':
                eval_df = eval_disease(disease, clinvar_df_pre_ls, disease_df[disease_df.Disease == disease],
                                       fout_stats, fout_eval, clin_labels, score_cols)
                eval_df_ls.append(eval_df)
    pd.concat(eval_df_ls).to_csv(args.out_df, index=False, sep='\t')

if __name__ == "__main__":
    desc = 'Eval global mpc cutoff'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('score_cols', 'clinvar', 'clinvar_single', 'clinvar_mult', 'clinvar_exp', 'denovo',
             'gene_dx', 'uc', 'other_disease', 'stats_out', 'eval_out', 'out_df')
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
