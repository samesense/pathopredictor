"""Handle neutral uniprot variants for FATHMM.
   I onlyget have codons here, so you are in the
   neutral dataset if you fall in the codon.
"""

up_ftp = 'ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds/UP000005640_9606_variants.bed'
rule dl_uniprot_humsavar:
    input: HTTP.remote("www.uniprot.org/docs/humsavar.txt", keep_local=True)
    output: DATA + 'raw/uniprot/humsavar'
    shell:  'mv {input} {output}'

rule dl_uniprot_translation:
    input: FTP.remote(up_ftp, keep_local=True)
    output: DATA + 'raw/uniprot/variants.hg38.bed'
    shell: 'mv {input} {output}'

rule add_pos_to_uniprot:
    input:  uniprot = DATA + 'raw/uniprot/humsavar',
            pos = DATA + 'raw/uniprot/variants.hg38.bed'
    output: o = DATA + 'interim/uniprot/humsavar.hg38.pos',
            miss = DATA + 'interim/uniprot/humsavar.miss'
    run:
        vars = {}
        read_vars = False
        with open(input.uniprot) as f:
            for line in f:
                if read_vars:
                    sp = line.strip().split()
                    if not sp:
                        break
                    if sp[4] == 'Polymorphism':
                        vars[ sp[2] ] = True
                elif 'gene' == line[:4]:
                    f.readline()
                    read_vars = True
        with open(input.pos) as f, open(output.o, 'w') as fout:
            for line in f:
                sp = line.strip().split('\t')
                var = sp[12]
                if var in vars:
                    # they give the whole codon in bed format
                    for pos in range(int(sp[1])+1, int(sp[2])+1):
                       ls = (sp[0], str(pos))
                       print('\t'.join(ls), file=fout)
                    vars[var] = True
        with open(output.miss, 'w') as fout:
            for var in vars:
                if not vars[var]:
                    print(var, file=fout)
