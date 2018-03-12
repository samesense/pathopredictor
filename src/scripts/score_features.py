"""Calculate average precision, auc
"""
import argparse
import pandas as pd
from sklearn import metrics

def mk_pred_col(features):
    if '-' in features:
        return features + '_probaPred'
    else:
        return features

def run_anova(features, df):
    pval_ls = []
    for feature in features.split('-'):
        pval_ls.append(anova(features, feature, df))
    return max(pval_ls)

def main(features, in_file, out_file):
    df = pd.read_csv(in_file, sep='\t')
    pred_col = mk_pred_col(features)
    fpr, tpr, _ = metrics.roc_curve(df['y'], df[pred_col], pos_label=1)
    auc = metrics.auc(fpr, tpr)
    avg_pr = metrics.average_precision_score(df['y'], df[pred_col])

    out_df = pd.DataFrame({'features':features,
                           'auc':auc, 'avg_pr':avg_pr,})
    out_df.to_csv(out_file, index=False, sep='\t')

if __name__ == "__main__":
    features, in_file, out_file = sys.argv[1:]
    main(features, in_file, out_file)
