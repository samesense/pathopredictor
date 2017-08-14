import pandas, sys
dat_file, vcf_file, out_file = sys.argv[1:]
df_pre = pandas.read_excel(dat_file)
genes = set(df_pre['Gene Symbol'].values)

with open(vcf_file) as f, open(out_file, 'w') as fout:
    for line in f:
        if line[0] == '#':
            print(line.strip(), file=fout)
        else:
            sp = line.strip().split('\t')
            effs = sp[-3].split('EFF=')[1].split(';')[0].split(',')
            for eff in effs:
                # downstream_gene_variant(MODIFIER||958|||WASH7P||NON_CODING|NR_024540.1||1)
#                print(eff)
                gene = eff.split('|')[-6]
                if gene in genes:
                    print(line.strip(), file=fout)
