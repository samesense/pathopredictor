"""Discard variants used to train classifiers used as features."""

rule combine_discard:
    input: hgmd = ,
           uniprot=
    output: DATA + 'processed/discard_variants'
