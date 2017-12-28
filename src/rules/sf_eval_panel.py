include: "const.py"

rule eval_panel_global:
    input:  DATA + 'interim/clinvar/clinvar.limit3.dat',
            DATA + 'interim/EPIv6.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/uc.eff.dbnsfp.anno.hHack.dat.limit.xls',
            DATA + 'interim/other/other.eff.dbnsfp.anno.hHack.dat.limit.xls'
    output: WORK + 'eval_panel.stats',
            WORK + 'eval_panel.eval'
    shell:  'python {SCRIPTS}score_panel_global_model.py {input} {output}'

rule limit_for_plot:
    input:  WORK + 'eval_panel.eval'
    output: WORK + 'global.eval_panel.eval.totWrong'
    shell:  "grep 'TotWrong\|disease' {input} | grep -v ssue | grep -v earing > {output}"

# ggplot(data=d) + geom_col(aes(y=var_count,x=score_type, fill=score_type)) + facet_grid(disease~., scale='free') + theme_bw() + ylab('Wrong Predictions') + theme(axis.text.x = element_text(angle=90)) + xlab('') + theme(legend.position="none")    


