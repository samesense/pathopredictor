"""Mk data dumps for manuscript."""

rule paper_training_data:
    input:  i = DATA + 'interim/full/panel.dat'
    output: o = DATA + 'processed/dryad/S1_trainingData_hg19.csv'
    run:
        cols = ['chrom', 'pos', 'ref', 'alt', 'class', 'Disease']
        pd.read_csv(input.i, sep='\t')[cols].to_csv(output.o, index=False, sep=',')

rule upload_supp:
    input:  DATA + 'processed/dryad/{table}.csv'
    output: DBox.remote('ahmad_predictor/{table}.csv')
    shell:  'cp {input} {output}'
