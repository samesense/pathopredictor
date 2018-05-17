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
