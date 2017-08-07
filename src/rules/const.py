import os, sys

SECRETS = '/home/evansj/me/.secrets/'
sys.path.append(SECRETS)
from pass_wd import *

DONE = TWILIO_PRE + "--data-urlencode 'Body=DONE' " + TWILIO_POST
FAIL = TWILIO_PRE + "--data-urlencode 'Body=FAIL' " + TWILIO_POST

p = os.getcwd()
if 'src' in p:
    PWD = p.split('src/')[0]
else:
    PWD = p + '/'

WORK = PWD + 'work/'
FILES = PWD + 'docs/'
SCRIPTS = PWD + 'src/scripts/'
DATA = PWD + 'data/'
LOG = PWD + 'log/'
CONFIG = PWD + 'configs/'
