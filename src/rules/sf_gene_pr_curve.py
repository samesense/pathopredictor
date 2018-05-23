"""Gene-wise PR curves for top genes after global model"""

def calc_pr_curve(rows, gene):
    pred_col = [x for x in rows.columns if 'probaPred' in x][0]
    precision, recall, thresholds = precision_recall_curve(rows['y'], rows[pred_col], pos_label=1)
    fpr, tpr, thresholds = roc_curve(rows['y'], rows[pred_col], pos_label=1)
    curve = {'fpr':fpr, 'tpr':tpr}
    s = pd.DataFrame(curve, columns=['pre', 'rec', 'fpr', 'tpr'])
    s.loc[:, 'gene'] = gene
    return s

rule calc_pr_curve:
    input:  i = WORK + 'roc_df_{dat}/{cols}'
    output: o = DATA + 'interim/gene_pr/{dat}.{cols}'
    run:
        df = pd.read_csv(input.i, sep='\t')
        ls = []
        for gene in set(df['gene']):
            dd = df[df.gene==gene]
            ls.append(calc_pr_curve(dd, gene))
        pd.concat(ls).to_csv(output.o, index=False, sep='\t')

rule plot_pr_curve:
    input: p = DATA + 'interim/gene_pr/panel.' + C_FEATS,
           c = DATA + 'interim/gene_pr/clinvar.' + C_FEATS
    output: o = DOCS + 'paper_plts/fig6_byGene.pdf'
    run:
        genes = ['KCNQ2', 'STXBP1',
                 'SCN2A', 'SCN5A', 'RAF1']
        panel = pd.read_csv(input.p, sep='\t')
        p_crit = panel.apply(lambda row: row['gene'] in genes and row['gene'] != 'RAF1', axis=1)
        clinvar = pd.read_csv(input.c, sep='\t')
        c_crit = clinvar.apply(lambda row: row['gene'] in genes, axis=1)
        panel.loc[:, 'dataset'] = 'Disease panel'
        clinvar.loc[:, 'dataset'] = 'Total ClinVar'
        pd.concat([panel[p_crit], clinvar[c_crit]]).to_csv(output.o + '.dat', index=False, sep='\t')
        R("""require(ggplot2)
             d = read.delim("{output}.dat", sep='\t', header=TRUE)
             p = ggplot(data=d) + geom_line(aes(x=fpr, y=tpr, colour=gene)) +
                 facet_grid(dataset~.) + theme_bw(base_size=18) +
                 xlab('False positive rate') + ylab('True positive rate')
            ggsave("{output}", p)
          """)

rule all_pr:
    input: DOCS + 'paper_plts/fig6_byGene.pdf'
