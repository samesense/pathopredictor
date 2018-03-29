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

def clinvar_stats(disease, clinvar_df_pre, disease_df, clinvar_label, gnomad_clinvar):
    """Remove disease panel testing data from clinvar"""
    key_cols = ['chrom', 'pos', 'ref', 'alt']
    test_keys = {':'.join([str(x) for x in v]):True for v in disease_df[key_cols].values}

    crit = clinvar_df_pre.apply(lambda row: not ':'.join([str(row[x]) for x in key_cols]) in test_keys, axis=1)
    clinvar_df = clinvar_df_pre[crit]

    disease_genes = set(disease_df['gene'])
    crit = clinvar_df.apply(lambda row: row['gene'] in disease_genes, axis=1)
    clinvar_df_limit_genes = clinvar_df[crit]

    # pad w/ benign from gnomad
    genes = set(clinvar_df_limit_genes['gene'])
    gnomad_dfs = []
    for gene in genes:
        crit = clinvar_df_limit_genes.apply(lambda row: row['gene'] == gene and row['y']==1, axis=1)
        path_count = len(clinvar_df_limit_genes[crit])

        crit = clinvar_df_limit_genes.apply(lambda row: row['gene'] == gene and row['y']==0, axis=1)
        benign_count = len(clinvar_df_limit_genes[crit])

        if benign_count < path_count:
            rows = len(gnomad_clinvar[gnomad_clinvar.gene==gene])
            n = min([rows, path_count-benign_count])
            print(gene, n, path_count, benign_count)
            if rows:
                df = gnomad_clinvar[gnomad_clinvar.gene==gene].sample(n)
                gnomad_dfs.append(df)
    if gnomad_dfs:
        l = pd.concat([clinvar_df_limit_genes] + gnomad_dfs)
    else:
        l = clinvar_df_limit_genes
    # clinvar_df not used
    return clinvar_df, l

def print_data_stats(disease, clinvar_df_pre_ls, disease_df, clin_labels, gnomad_clinvar):
    clin_dat = list(map(lambda x: clinvar_stats(disease, x[0], disease_df, x[1], gnomad_clinvar),
                        list(zip(clinvar_df_pre_ls, clin_labels))))

    return [x[0] for x in clin_dat], [x[1] for x in clin_dat]

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

def eval_disease(disease, clinvar_df_pre_ls, disease_df, clin_labels, cols, gnomad_clinvar, gnomad_panel):
    _, clinvar_df_limit_genes_ls = print_data_stats(disease, clinvar_df_pre_ls, disease_df, clin_labels, gnomad_clinvar)
    # keep track of clinvar results# keep track of clinvar results
    for label, df in zip(clin_labels, clinvar_df_limit_genes_ls):
        df['Disease'] = disease + ':' + label

    # one gene at a time
    clinvar_acc = []
    genes = set(disease_df['gene'])

    print(disease, 'one', disease_df.columns.values)
    for test_gene in genes:
        # pad disease w/ benign
        df = disease_df[disease_df.gene==test_gene]
        print(df.columns.values)
        crit = df.apply(lambda row: row['y']==1, axis=1)
        path_count = len(df[crit])

        crit = df.apply(lambda row: row['y']==0, axis=1)
        benign_count = len(df[crit])

        gnomad_dfs = []
        if benign_count < path_count:
            rows = len(gnomad_panel[gnomad_panel.gene==test_gene])
            n = min([rows, path_count-benign_count])
            if rows:
                dfg = gnomad_panel[gnomad_panel.gene==test_gene].sample(n)
                gnomad_dfs.append(dfg)
        if path_count < 5:
            continue

        if gnomad_dfs:
            sub_train_df = pd.concat([df] + gnomad_dfs)
        else:
            sub_train_df = df

        crit = sub_train_df.apply(lambda row: row['y']==0, axis=1)
        benign_count = len(sub_train_df[crit])
        if benign_count < 5:
            continue

        tree_clf_sub = tree.DecisionTreeClassifier( max_depth=len(cols) )
        X, y = sub_train_df[cols], sub_train_df['y']
        sub_train_df.to_csv('shit.out', index=False, sep='\t')
        tree_clf_sub.fit(X, y)

        lm =  linear_model.LogisticRegression(fit_intercept=True)
        lm.fit(X, y)
        pred_col = '-'.join(cols) + '_pred_lm'
        for clin_label, clinvar_df_full in zip(clin_labels, clinvar_df_limit_genes_ls):
            # limit to gene
            clinvar_df = clinvar_df_full[clinvar_df_full.gene==test_gene]
            X = clinvar_df[cols]
            if len(clinvar_df):
                clinvar_preds = lm.predict(X) #tree_clf_sub.predict(X)
                lm_preds = lm.predict(X)
                clinvar_df.loc[:, pred_col] = lm_preds
                lm_proba = [x[1] for x in lm.predict_proba(X) ]
                clinvar_df.loc[:, '-'.join(cols) + '_probaPred'] = lm_proba
                #clinvar_df['mpc_pred'] = clinvar_preds
                clinvar_df.loc[:, 'PredictionStatusMPC'] = clinvar_df.apply(lambda row: eval_pred(row, pred_col), axis=1)
                clinvar_df.loc[:, 'PredictionStatusBaseline'] = clinvar_df.apply(lambda x: eval_mpc_raw(x, cols), axis=1)
                clinvar_acc.append(clinvar_df)

    if clinvar_acc:
        clinvar_result_df = pd.concat(clinvar_acc)
    else:
        clinvar_result_df = pd.DataFrame()
    return clinvar_result_df

def main(args):
    score_cols = args.score_cols.split('-')
    FOCUS_GENES = ('SCN1A','SCN2A','KCNQ2', 'KCNQ3', 'CDKL5',
                   'PCDH19', 'SCN1B', 'SCN8A', 'SLC2A1',
                   'SPTAN1', 'STXBP1', 'TSC1')

    # load clinvar benign padding
    gnomad_clinvar = pd.read_csv(args.gnomad_clinvar, sep='\t')
    gnomad_clinvar.loc[:, 'is_domain'] = gnomad_clinvar.apply(lambda row: 0 if 'none' in row else 1, axis=1)
    gnomad_clinvar['y'] = 0

    # load panel benign padding
    gnomad_panel= pd.read_csv(args.gnomad_clinvar, sep='\t')
    gnomad_panel.loc[:, 'is_domain'] = gnomad_panel.apply(lambda row: 0 if 'none' in row else 1, axis=1)
    gnomad_panel['y'] = 0

    # load clinvar
    clin_labels = ('tot', 'single', 'mult', 'exp',)
    clinvars = (args.clinvar, args.clinvar_single, args.clinvar_mult, args.clinvar_exp)
    clinvar_df_pre_ls = list(map(lambda x: pd.read_csv(x, sep='\t').rename(columns={'clin_class':'y'}), clinvars))

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

    eval_df_clinvar_ls = []
    for disease in diseases:
        if str(disease) != 'nan' and not 'Connective' in disease:
            eval_clinvar_df = eval_disease(disease, clinvar_df_pre_ls, disease_df[disease_df.Disease == disease],
                                           clin_labels, score_cols, gnomad_clinvar, gnomad_panel)
            eval_df_clinvar_ls.append(eval_clinvar_df)
    pd.concat(eval_df_clinvar_ls).to_csv(args.out_df_clinvar_eval, index=False, sep='\t')

if __name__ == "__main__":
    desc = 'Eval combinations of features on panel and clinvar'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('score_cols', 'clinvar', 'clinvar_single', 'clinvar_mult', 'clinvar_exp',
             'gene_dx', 'uc', 'other_disease', 'gnomad_panel', 'gnomad_clinvar', 'out_df_clinvar_eval')
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
