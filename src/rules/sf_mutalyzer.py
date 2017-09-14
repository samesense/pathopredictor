"""Mk input for mutalyzer"""
from const import *

rule panel_one:
    input:  DATA + 'raw/EPIv6.xlsx'
    output: DATA + 'interim/mutalyzer_input/panel_one'
    shell:  'python {SCRIPTS}mk_dat.py {input} {output}'

rule panel_two:
    input:  DATA + 'raw/EpilepsyVariantDataForAhmadClean_090517.xlsx'
    output: DATA + 'interim/mutalyzer_input/panel_two'
    shell:  'python {SCRIPTS}mk_dat_panel_two.py {input} {output}'
