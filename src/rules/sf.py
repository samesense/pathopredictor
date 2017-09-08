"""Main snakefile"""

include: "const.py"
include: "sf_ann.py"
include: "sf_clinvar.py"

rule all:
    input: DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.xls', DATA + 'interim/clinvar/clinvar.dat'
