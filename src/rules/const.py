import os, sys
# sys.path.append('/home/evansj/me/projects/me/tool_dirs/')
# from tools import *

# SECRETS = '/home/evansj/me/.secrets/'
# sys.path.append(SECRETS)
# from pass_wd import *

# DONE = TWILIO_PRE + "--data-urlencode 'Body=DONE' " + TWILIO_POST
# FAIL = TWILIO_PRE + "--data-urlencode 'Body=FAIL' " + TWILIO_POST

p = os.getcwd()
if 'src' in p:
    PWD = p.split('src/')[0]
else:
    PWD = p + '/'

TMP = PWD + 'tmp/'
WORK = PWD + 'work/'
WORK2 = '/mnt/isilon/cbmi/variome/perry/projects/sarmadi/mahdi_epi/work/'
DOCS = PWD + 'docs/'
FILES = PWD + 'docs/'
SCRIPTS = PWD + 'src/scripts/'
SCRIPTS2 = '/mnt/isilon/cbmi/variome/perry/projects/sarmadi/mahdi_epi/src/scripts/'
DATA = PWD + 'data/'
LOG = PWD + 'log/'
CONFIG = PWD + 'configs/'
PLOTS = PWD + 'docs/plots/'

EXAC_DIR = '/home/evansj/me/projects/diskin/target_exac_setup/data/'
EXAC_PED = '/home/evansj/me/projects/diskin/target_exac_setup/files/JUNK_PED.ped'
VCFANNO_LUA_FILE = '/home/evansj/me/projects/me/vcfanno_lua/scripts/target.lua'

ISILON = '/mnt/isilon/cbmi/variome/perry/'
TABIX = ISILON + 'bin/tabix'
BGZ = '/nas/is1/bin/bgzip'
VT = ISILON + 'condas/miniconda3/envs/target_capture/bin/vt'
PY27_T = '~/me/condas/miniconda3/envs/testPy27/bin/python'
PY27 = '~/me/condas/miniconda3/envs/py27/bin/python'
JAVA = '/nas/is1/bin/java'
EFF = '/home/evansj/me/tools/snpEff/snpEff.jar'
SIFT = '/home/evansj/me/tools/snpEff/SnpSift.jar'
SIFT_DBNSFP = '/home/evansj/me/tools/snpEff/data/dbNSFP/dbNSFP2.4.txt.gz'
EFF_CONFIG = '/home/evansj/me/tools/snpEff/snpEff.config'
CLINVAR = '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/clinvar_20160203_noUnicode.tidy.vcf.gz'
HEADER_HCKR =  '/nas/is1/perry/projects/me/vcfHeaderHckr/vcfHeadrHckr.py'
VCFANNO = '/mnt/isilon/cbmi/variome/bin/vcfanno/0.0.11/bin/vcfanno'
GEMINI_ANNO = '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data'

HG19_FA = DATA + 'raw/hg19.fa'
