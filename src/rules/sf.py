"""Main snakefile"""

include: "const.py"

rule all:
    output: LOG + 'DONE'
    shell:  'touch {output}'
