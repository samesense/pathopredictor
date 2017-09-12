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

WORK = PWD + 'work/'
DOCS = PWD + 'docs/'
FILES = PWD + 'docs/'
SCRIPTS = PWD + 'src/scripts/'
DATA = PWD + 'data/'
LOG = PWD + 'log/'
CONFIG = PWD + 'configs/'
PLOTS = PWD + 'docs/plots/'

EXAC_DIR = '/home/evansj/me/projects/diskin/target_exac_setup/data/'
EXAC_PED = '/home/evansj/me/projects/diskin/target_exac_setup/files/JUNK_PED.ped'
VCFANNO_LUA_FILE = '/home/evansj/me/projects/me/vcfanno_lua/scripts/target.lua'

ISILON = '/mnt/isilon/cbmi/variome/perry/'
TABIX = ISILON + 'bin/tabix'
BGZ = ISILON + 'bin/bgzip'
VT = ISILON + 'condas/miniconda3/envs/target_capture/bin/vt'
PY27_T = '~/me/condas/miniconda3/envs/testPy27/bin/python'

CLINVAR = '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/clinvar_20160203_noUnicode.tidy.vcf.gz'
