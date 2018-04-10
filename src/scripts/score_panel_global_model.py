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
    """Remove disease panel testing data from clinvar"""
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
        #print(X_test)
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

def eval_disease(disease, clinvar_df, disease_df, cols):
    clin_labels = ('clinvar_tot', 'clinvar_single')
    # one gene at a time
    acc_df_ls = []
    clinvar_acc = []
    genes = set(disease_df['gene'])

    for test_gene in genes:
        sub_train_df = disease_df[disease_df.gene != test_gene]
        tree_clf_sub = tree.DecisionTreeClassifier( max_depth=len(cols) )
        X, y = sub_train_df[cols], sub_train_df['y']
        tree_clf_sub.fit(X, y)

        lm =  linear_model.LogisticRegression(fit_intercept=True)
        lm.fit(X, y)

        test_df = disease_df[disease_df.gene == test_gene]
        X_test = test_df[cols]
        preds = tree_clf_sub.predict(X_test)
        lm_preds = lm.predict(X_test)
        pred_col = '-'.join(cols) + '_pred_lm'
        test_df.loc[:, pred_col] = lm_preds
        lm_proba = [x[1] for x in lm.predict_proba(X_test) ]
        test_df.loc[:, '-'.join(cols) + '_probaPred'] = lm_proba
        test_df.insert(2, 'PredictionStatus', test_df.apply(lambda row: eval_pred(row, pred_col), axis=1) )

        acc_df_ls.append(test_df)

        for clin_label in clin_labels:
            # limit to gene
            clinvar_df = clinvar_df[(clinvar_df.gene==test_gene) & (clinvar_df.clinvar_subset==clin_label)]
            clinvar_df['Disease'] = disease
            X = clinvar_df[cols]
            if len(clinvar_df):
                clinvar_preds = lm.predict(X)
                lm_preds = lm.predict(X)
                clinvar_df.loc[:, pred_col] = lm_preds
                lm_proba = [x[1] for x in lm.predict_proba(X) ]
                clinvar_df.loc[:, '-'.join(cols) + '_probaPred'] = lm_proba
                clinvar_df.loc[:, 'PredictionStatus'] = clinvar_df.apply(lambda row: eval_pred(row, pred_col), axis=1)
                clinvar_acc.append(clinvar_df)

    test_df = pd.concat(acc_df_ls)
    clinvar_result_df = pd.concat(clinvar_acc)
    return test_df, clinvar_result_df

def main(args):
    score_cols = args.score_cols.split('-')
    FOCUS_GENES = ('SCN1A','SCN2A','KCNQ2', 'KCNQ3', 'CDKL5',
                   'PCDH19', 'SCN1B', 'SCN8A', 'SLC2A1',
                   'SPTAN1', 'STXBP1', 'TSC1')

    clinvar_df = pd.read_csv(args.clinvar, sep='\t').rename(columns={'Disease':'clinvar_subset'})

    # load panels
    panel_df = pd.read_csv(args.panel, sep='\t')
    crit = panel_df.apply(lambda row: row['gene'] in FOCUS_GENES, axis=1)
    disease_genedx_limitGene_df= panel_df[crit]
    disease_genedx_limitGene_df['Disease'] = 'genedx-epi-limitGene'

    disease_df = pd.concat([panel_df, disease_genedx_limitGene_df])
    diseases = set(disease_df['Disease'])

    eval_df_ls, eval_df_clinvar_ls = [], []
    for disease in diseases:
        eval_panel_df, eval_clinvar_df = eval_disease(disease, clinvar_df, disease_df[disease_df.Disease == disease],
                                                      score_cols)
        eval_df_ls.append(eval_panel_df)
        eval_df_clinvar_ls.append(eval_clinvar_df)
    pd.concat(eval_df_ls).to_csv(args.out_df, index=False, sep='\t')
    pd.concat(eval_df_clinvar_ls).to_csv(args.out_df_clinvar_eval, index=False, sep='\t')

if __name__ == "__main__":
    desc = 'Eval combinations of features on panel and clinvar'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('score_cols', 'clinvar', 'panel',
             'out_df', 'out_df_clinvar_eval')
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
