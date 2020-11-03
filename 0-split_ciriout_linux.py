#!usr/bin/python3
from pathlib import Path
import os, sys, re,shutil,multiprocessing,time
maindir=str(Path(__file__).parent)
os.chdir(maindir)
split_num=88

if os.path.isdir(maindir+"/tmp"):
    shutil.rmtree(maindir+"/tmp")
time.sleep(2)
for num in range(split_num):
    tmpname=maindir+"/tmp/split{num}".format(num=num)
    os.makedirs(tmpname)
time.sleep(1)
with open(maindir+'/CIRI.result') as infile:
    files = [open(maindir+'/tmp/split{num}/split{num}.result'.format(num=num), 'w') for num in range(0,split_num)]
    for linenum, linecont in enumerate(infile):
        files[linenum % split_num].write(linecont)
    for spf in files:
        spf.close()
time.sleep(2)
