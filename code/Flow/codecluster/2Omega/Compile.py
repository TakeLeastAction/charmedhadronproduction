#!/usr/bin/env python
#coding=gbk
import glob
import os
import time
import datetime
import subprocess
import re
import string
import math
from shutil import copy
import smtplib
from email.mime.text import MIMEText
from email.header import Header
from scipy.optimize import fmin_cobyla     # Optimization method


cmd='ifort -O1 col_xn.for -o col_xn'
os.system(cmd)
p1=open('jobNum.txt')
p2=p1.read()
jobNumL=[]
jobNumL.append(string.atof(p2.strip(' \n')))
jobNum=int(jobNumL[0])
for i in range(1,jobNum+1):
    if i<10:
        kdr='00'+str(i)
        os.mkdir(kdr)
    elif i<100:
        kdr='0'+str(i)
        os.mkdir(kdr)
    else:
        kdr=str(i)
        os.mkdir(kdr)
        
    copy('input.txt',kdr)
    copy('jobNum.txt',kdr)
    copy('PT.txt',kdr)
    copy('col_xn',kdr)
    os.chdir(kdr)
    ss=open('circle.txt','w')
    print >>ss,i
    ss.close()  
    os.chdir('..')
        
        
        
        
