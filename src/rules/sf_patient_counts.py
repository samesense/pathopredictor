"""Counts of patients on panels"""

rule epi_counts:
    input:  i = DATA + 'raw/EPIv6.xlsx'
    output: o = DATA + 'interim/patient_counts/epi'
    run:
        df = pandas.read_excel(input.i)
        df.loc[:, 'final_count'] = df.apply(lambda row: row['Pos Fam Cnt'] + row['Neg Fam Cnt'], axis=1)
        max_count = max(df['final_count'])
        with open(output.o, 'w') as fout:
            print('Epilesy\t' + str(max_count), file=fout)

rule other_panel_counts:
    input:  i = DATA + 'raw/Clinicalvariants_LMM.xlsx',
            gene_info = DATA + 'raw/gene_disease.xlsx'
    output: o = DATA + 'interim/patient_counts/other'
    run:
        df = pd.read_excel(input.i).rename(columns={'Gene Symbol':'gene'})
        disease_df = pd.read_excel(input.gene_info, skiprows=[0,1,2]).rename(columns={'Gene':'gene'})
        df = pd.merge(disease_df, df, on='gene', how='left')
        df.groupby('Disease')[['# of probands']].max().reset_index().to_csv(output.o, index=False, sep='\t')
        
