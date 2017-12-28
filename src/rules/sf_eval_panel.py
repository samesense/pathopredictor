include: "const.py"

rule eval_panel_global:
    input:  DATA + 'interim/clinvar/clinvar.limit3.dat',
            DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/other/other.eff.dbnsfp.anno.hHack.dat.limit.xls'
    output: WORK + 'eval_panel.stats',
            WORK + 'eval_panel.eval'
    shell:  'python {SCRIPTS}score_panel_global_model.py {input} {output}'            
