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

totalJob=300
sptime=1
sleeptime=0.5

def judgeDir():
    dirfile="1"
    while True:
        if os.path.isdir("./"+dirfile) :
            return len(dirfile)
        else :
            dirfile="0"+dirfile
            
def changeint(num,length):    # num=1,length=3; "001"
    string=str(num)
    return "0"*(length-len(string))+string

def jobNum():
    jobNum=0
    f=subprocess.Popen(("condor_q"),stdout=subprocess.PIPE).stdout.readlines()
    
    for i in f:
        if re.search("sunkj.* ",i):
            jobNum+=1
            
  
          
    print 'the number of jobs'     
    print jobNum        
    return jobNum

def getDirList( dirpath ):
    s1=os.getcwd()
    os.chdir(dirpath)
    s=os.getcwd()
    L = os.listdir(s)
    p=str(s)
    b = [ p +'/'+ x    for x in L if os.path.isdir( p +'/'+ x ) ]
    os.chdir(s1)
    return b


def submitjob(dir):
   
    os.chdir(dir)
    cmd='./'+'/makejob'    
    os.system(cmd)
    s=os.getcwd()
    L=os.listdir(s)
    print L
    p=str(s)
    b = [ x   for x in L if os.path.isdir( p +'/'+ x ) ]
    print b
    copy('input.txt',b[0])
    copy('circle.txt',b[0])
    os.system('./'+'makejob-condor')
    os.chdir('..')
    
def getJobList(dirpath,jobname):
    s=os.getcwd()
    os.chdir(dirpath)
    ls = [os.path.join(root, name)
    for root, dirs, files in os.walk(os.getcwd())
        for name in files
            if name.endswith(jobname)]
    os.chdir(s)    
    return ls            
        

def sendMail(res,runtime):
    sender = 'dcbasunsiabc@163.com' 
    receiver = 'sunkaijia@sjtu.edu.cn'
    subject = 'working report from server'
    smtphost = 'smtp.163.com'
    username = 'dcbasunsiabc'
    password = 'abc667667'
    message1 = "The result is the "+str(res)
    message2 = " Lets cheers "
    message3 = "The running time is "+str(runtime)
    message = message1+message2+message3
    msg = MIMEText(message,'text','utf-8')
    msg['Subject'] = Header(subject,'utf-8')
    smtp = smtplib.SMTP()
    smtp.connect('smtp.163.com')
    smtp.login(username, password)
    smtp.sendmail(sender, receiver, msg.as_string())
    smtp.quit()
           
  
def throw_job(dirpath):        
    for i in dirpath:
        print os.getcwd()
        print i
        copy('input.txt',i)
        os.chdir(i)
        cmd0='ls|grep [0-9]|xargs rm -rf'
        cmd1='chmod'+' +x'+' play.sh'  
        cmd2='./play.sh'
        os.system(cmd0)
        os.system(cmd1)
        os.system(cmd2)
        time.sleep(sptime)
        os.chdir('..')   

def status_job(dirpath,jobname):
    s=''
    for i in range(0,len(dirpath)):
        dirlis=getDirList( dirpath[i] )
        joblis=getJobList( dirpath[i],jobname[i] )
        if len(dirlis)>len(joblis):
            s= s + dirpath[i]
            
    return s        
    

def getData(dirpath,jobname):
    dirlis=getJobList( dirpath,jobname )
    dirlis=sorted(dirlis)
    spectrum=[]
    for i in dirlis:
        p1=open(i)
        p2=p1.read()
        spectrum.append(p2.strip(' \n'))
        
    return spectrum

def saveData(jobName,jobData):
    ss=open(jobName,'w')
    for i in range(0,len(jobData)):
        print >>ss,jobData[i]
    ss.close()     

def objFun():
    # 将参数复制到文件夹(们) b 
    s=os.getcwd()
    L=os.listdir(s)
    print L
    p=str(s)
    b = [ p+'/1phi',p+'/2Omega' ]
    basic_dir=b
    Particle=['phi.txt','Omega.txt']
    print b
    # 进入到每个文件夹去，提交 fortran 程序 得到 
    det=0
    while det==1:
        throw_job(b) 
        jobnum=jobNum()   
        print jobnum
        while jobnum >= 1:
            time.sleep(sleeptime)
            jobnum=jobNum()  
  
        det=0          
       # pathdir=status_job(b,Particle) 
       # b=pathdir
       # if (len(b)==0):
       #     det=0
  
    # phi
    Data_phi=getData(basic_dir[0],Particle[0])
    saveData("./result_phi.dat",Data_phi)   
    # Omega
    Data_Omega=getData(basic_dir[1],Particle[1])
    saveData("./result_Omega.dat",Data_Omega) 
    # proton
    #Data_p=getData(basic_dir[2],Particle[2])
    #saveData("./result_p.txt",Data_p) 
    # Lambda
    #Data_lambda=getData(basic_dir[3],Particle[3])
    #saveData("./result_lambda.txt",Data_lambda) 
    # Xi
    #Data_xi=getData(basic_dir[4],Particle[4])
    #saveData("./result_xi.txt",Data_xi)     
    # Occc
    #Data_Occc=getData(basic_dir[5],Particle[5])
    #saveData("./result_Occc.txt",Data_Occc)
    
objFun()  
