#!/usr/bin/python
#####################################################
#
# Baldwin standard library
#
# 21st August 2012
####################################################


#use countmaxthread in the path to return the max number of 
#threads, using openMp parallelisation
def GetMaxThreads():
    import subprocess as sp,string,sys
    p1=sp.Popen(["CountMaxThread"],stdout = sp.PIPE)
    nthreads=int(string.split(string.split(p1.stdout.read(),'\n')[0])[1])
    sys.stdout.write('Detected maximum of %i threads on system\n' % nthreads)
    return nthreads
    
#determine if a path exists. If not, create it.
def PathExists(pathList):
    import os,sys
    for i in range(len(pathList)):
        if(os.path.exists(pathList[i])==0):
            sys.stdout.write('Creating %s\n' % pathList[i])
            os.mkdir(pathList[i])
        else:
            #sys.stdout.write('Path already exists: %s\n' % pathList[i])
            pass

#read a file and return a space delimited array
def readfile(infile):
    import string
    peak=[]
    peakfile=open(infile,'r')
    for line in peakfile.readlines():
        linetosave=string.split(line)
        peak.append(linetosave)
    peakfile.close()
    return peak
