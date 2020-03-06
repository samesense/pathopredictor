"""Annotate epilepsy panel vcfs"""

rule snpeff_general:
    input:  vcf = DATA + 'interim/epi/merged_all.vcf.gz',
    output: DATA + 'interim/{lab_dir,epi|other|user_preds}/{lab}.eff.vcf'
    shell:  """java -Xmx32g -Xms16g -jar {EFF} eff -dataDir {DOCKER_DATA}snpeff/data/ \
               -strict -noStats GRCh37.75 -c {EFF_CONFIG} \
               {input} > {output}"""

rule annotateDbnsfp_general:
    input:  vcf = DATA + 'interim/{dir}/{lab}.eff.vcf',
            db = DATA + 'raw/snpsift/dbNSFP.txt.gz'
    output: vcf = DATA + 'interim/{dir}/{lab}.eff.dbnsfp.tmp.vcf',
    shell:  """java -Xmx32g -Xms16g -jar {SIFT} dbnsfp \
               -db {input.db} -f {DBNSFP_FIELDS} {input.vcf} > {output}"""

rule fix_dbnsfp_general:
    input:  i = DATA + 'interim/{dir}/{lab}.eff.dbnsfp.tmp.vcf'
    output: o = DATA + 'interim/{dir}/{lab}.eff.dbnsfp.vcf'
    run:
        lines = False
        with open(input.i) as f, open(output.o, 'w') as fout:
            for line in f:
                if lines:
                    if line[0] == '#':
                        break
                print(line.strip(), file=fout)
                if line[0] != '#':
                    lines = True
# fix pfam
# /mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/hg19.pfam.ucscgenes.enum.bed.gz
# ann fixed pfam
# parse genes
# sf_predict and base pipeline need to split here to handle diff annotations
rule vcfanno_general:
    input:   vcf = DATA + 'interim/{dir}/{lab}.eff.dbnsfp.vcf',
             conf = CONFIG + 'vcfanno.predict.conf',
             lua = VCFANNO_LUA_FILE
    output:  DATA + 'interim/{dir}/{lab}.eff.dbnsfp.anno.vcf'
    threads: 10
    singularity: 'docker://szilvajuhos/sarek-vcfanno'
    shell:   """vcfanno -p {threads} -base-path {GEMINI_ANNO} -lua {input.lua} \
                {input.conf} {input.vcf} > {output}"""

#neg fam counts ppl
# pos fam counts ppl
# need to convert hom to a count of two
rule parse_vcf_general:
   input:  i = DATA + 'interim/{dir}/{lab}.eff.dbnsfp.anno.vcf'
   output: o = DATA + 'interim/{dir,user_preds|man|epi|other}/{lab}.eff.dbnsfp.anno.dat.xls',
   run:
       with open(input.i) as f, open(output.o, 'w') as fout:
           fields = ['chrom', 'pos', 'ref', 'alt',
                     'clin_class', 'pfam', 'eff', 'gene',
                     'esp_af_max', 'revel', 'mpc', 'mtr',
                     'ccr', 'fathmm', 'vest', 'missense_badness', 'missense_depletion']
           print('\t'.join(fields), file=fout)
           for line in f:
               if not line[0] == '#':
                   data = parse_vcf_data(line)
                   if data:
                       ls = [data[x] for x in fields]
                       print('\t'.join(ls), file=fout)

rule mahdi_ann:
    input:
        DATA + 'interim/epi/merged_all.eff.dbnsfp.anno.dat.xls'
        #DATA + 'interim/epi/merged_all.eff.dbnsfp.vcf'
