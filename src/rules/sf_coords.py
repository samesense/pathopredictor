"""Convert difficult transcripts"""

rule get_ncbi_seq:
    output: DATA + 'raw/fa/{nm}.fa'
    shell:  'python {SCRIPTS}get_ncbi_seq.py {wildcards.nm} {output}'

rule blat:
    input:  DATA + 'raw/fa/{nm}.fa'
    output: DATA + 'interim/blat/{nm}__blat'
    shell:  'blat {HG19_FA} {input} -out=pslx {output}'

rule coord_hash:
    input:  DATA + 'interim/blat/{nm}__blat',
            DATA + 'raw/fa/{nm}.fa'
    output: DATA + 'interim/blat_coord_hash/{nm}'
    shell:  'python {SCRIPTS}mk_coord_hash.py {input} {output}'

rule test_blat:
    input: expand( DATA + 'interim/blat_coord_hash/{nm}', nm = ('NM_001182.4', 'NM_000816.2', 'NM_000816.3') )
