"""Convert difficult transcripts.
   Run after mutalyzer.
"""

rule get_ncbi_seq:
    output: DATA + 'raw/fa/{nm}.cds.fa',
            DATA + 'raw/fa/{nm}.full.fa'
    shell:  'python {SCRIPTS}get_ncbi_seq.py {wildcards.nm} {output}'

rule blat:
    input:  DATA + 'raw/fa/{nm}.fa'
    output: DATA + 'interim/blat/{nm}__blat'
    shell:  'blat {HG19_FA} {input} -out=pslx {output}'

rule coord_hash:
    input:  DATA + 'interim/blat/{nm}.cds__blat',
            DATA + 'interim/blat/{nm}.full__blat',
            DATA + 'raw/fa/{nm}.cds.fa'
    output: DATA + 'interim/blat_coord_hash/{nm}'
    shell:  'python {SCRIPTS}mk_coord_hash.py {input} {output}'

def load_missing_transcripts(mutalyzer_file):
    """Find transcripts that failed mutalyzer"""
    nms = set()
    with open(mutalyzer_file) as f:
        f.readline()
        for line in f:
            sp = line.strip().split('\t')
            if len(sp) == 1:
                nms.add( sp[0].split(':')[0] )
    return nms

TS = load_missing_transcripts(DATA + 'raw/mutalyzer.panel_two.fix')

rule test_blat:
    input: expand( DATA + 'interim/blat_coord_hash/{nm}', nm = TS )
