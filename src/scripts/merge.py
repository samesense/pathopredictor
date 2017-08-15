"""Merge exac and fg"""
import pandas, math, argparse, sys
from collections import defaultdict

var_type = sys.argv[1]
fg_file = sys.argv[2] #'../data/interim/EPIv6.eff.dbnsfp.anno.hHack.splitPfam.dat'
bg_file = sys.argv[3] #'../data/interim/r1_no_tcga/exac.tidy.eff.dbnsfp.gt.anno.hHack.splitPfam.dat'
output = sys.argv[4]

df_exac_pre = pandas.read_csv(bg_file, sep='\t')
max_exac_an = max(df_exac_pre['an'].values)

df_fg_pre = pandas.read_csv(fg_file, sep='\t')
max_fg = max([p+2*n for p,n in df_fg_pre[['pos_fam', 'neg_fam']].values])
df_exac = df_exac_pre[df_exac_pre['af_1kg_all'] < .01]
df_fg = df_fg_pre[df_fg_pre['af_1kg_all'] < .01]

# don't forget to add non-var counts
var = var_type #'missense_variant'
cols = ['gene', 'pfam', 'chrom', 'pos']
pfam_to_exac_pos = defaultdict(dict)
for gene, pfam, chrom, pos in list(df_exac[df_exac.eff==var][cols].values):
    k = gene + 'xx' + pfam
    pfam_to_exac_pos[k][str(chrom) + ':' + str(pos)] = True
pfam_to_fg_pos = defaultdict(dict)
for gene, pfam, chrom, pos in list(df_fg[df_fg.eff==var][cols].values):
    k = gene + 'xx' + pfam
    pfam_to_fg_pos[k][str(chrom) + ':' + str(pos)] = True

missing_exac_count = []
for pfam in pfam_to_fg_pos:
    if not pfam in pfam_to_exac_pos:
        g,p = pfam.split('xx')
        ls = [g, p, len(pfam_to_fg_pos[pfam])]
        missing_exac_count.append(ls)
    else:
        miss_len = len( set(pfam_to_fg_pos[pfam]) - set(pfam_to_exac_pos[pfam]) )
        g,p = pfam.split('xx')
        ls = [g, p, miss_len]
        missing_exac_count.append(ls)
missing_exac_df = pandas.DataFrame({'gene':[x[0] for x in missing_exac_count],
                                    'pfam':[x[1] for x in missing_exac_count],
                                    'miss_exac':[x[2] for x in missing_exac_count]})        
missing_fg_count = []
for pfam in pfam_to_exac_pos:
    if not pfam in pfam_to_fg_pos:
        g,p = pfam.split('xx')
        ls = [g, p, len(pfam_to_exac_pos[pfam])]
        missing_fg_count.append(ls)
    else:
        miss_len = len( set(pfam_to_exac_pos[pfam]) - set(pfam_to_fg_pos[pfam]) )
        g,p = pfam.split('xx')
        ls = [g, p, miss_len]
        missing_fg_count.append(ls)
missing_fg_df = pandas.DataFrame({'pfam':[x[1] for x in missing_fg_count],
                                  'gene':[x[0] for x in missing_fg_count],
                                  'miss_fg':[x[2] for x in missing_fg_count]})

missing_df = pandas.merge(missing_fg_df, missing_exac_df, on=('gene', 'pfam'), how='outer').fillna(0)

#fg_pos = set([chrom + ':' + pos for chrom, pos in df_fg[df_fg.eff==var][cols].values])
#missing_exac_count = len(fg_pos) - len(exac_pos)

cols = ['ac', 'an']
g_exac = df_exac[df_exac.eff==var].groupby(['gene', 'pfam'])[cols].sum().reset_index()
#g_exac.head()

fg_cols = ['pos_fam', 'neg_fam']
g_fg = df_fg[df_fg.eff==var].groupby(['gene', 'pfam'])[fg_cols].sum().reset_index()
#g_fg.head()

m_pre = pandas.merge(g_fg, g_exac, on=('gene', 'pfam'), how='outer').fillna(0)
m = pandas.merge(m_pre, missing_df, on=('gene', 'pfam'), how='left')
print(max_exac_an, max_fg)
#m.head()

def calc_fg_neg(row, max_fg):
    return 2*row['neg_fam'] + max_fg*row['miss_fg'] + row['pos_fam']

def calc_bg_neg(row, max_exac):
    return row['an'] + max_exac*row['miss_exac']

m.loc[:, 'fg_tot'] = m.apply(lambda row: calc_fg_neg(row, max_fg), axis=1)
m.loc[:, 'bg_tot'] = m.apply(lambda row: calc_bg_neg(row, max_exac_an), axis=1)
m['fg_frac'] = (1+m['pos_fam']) / m['fg_tot']
m['fg_frac_log'] = m.apply(lambda row: math.log(row['fg_frac'], 2), axis=1)
m['bg_frac'] = (1+m['ac']) / m['bg_tot']
m['bg_frac_log'] = m.apply(lambda row: math.log(row['bg_frac'], 2), axis=1)
m['fg_other'] = m['fg_tot'] - m['pos_fam']
m['bg_other'] = m['bg_tot'] - m['ac']
m.to_csv(output, index=False, sep='\t')


