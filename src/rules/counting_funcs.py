import pandas as pd
def calc_tot_vars(rows):
    cols = ['tot_vars', 'CorrectPath', 'WrongPath', 'CorrectBenign', 'WrongBenign']
    tot_preds = sum(rows[rows.eval_type != 'TotWrong']['var_count'])
    counts = {row['eval_type']:int(row['var_count'])
              for _, row in rows.iterrows()
              if row['eval_type'] != 'TotWrong'}
    counts['tot_vars'] = int(tot_preds)
    cols_to_add = set(cols) - set(counts.keys())
    for c in cols_to_add:
        counts[c] = 0
    s = pd.Series(counts, index=cols)
    return s


