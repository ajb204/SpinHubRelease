import math, sys, os, string, random

######################################################
#statistical functions
######################################################
def covar(x,y):
    prod = 0
    n = len(x)
    for i in range(n):
        prod += x[i]*y[i]
    x_m = mean(x)
    y_m = mean(y)
    return prod/n - x_m*y_m
def mean(x):
    mean = 0
    n = len(x)
    for i in range(n):
        mean += x[i]
    return mean/n

def var(x):
    var = 0
    n = len(x)
    for i in range(n):
        var += (x[i]-mean(x))**2
    return var/n

def stdev(x):
    var = 0
    n = len(x)
    for i in range(n):
        var += (x[i]-mean(x))**2
    return (var/n)**0.5



def rmsd_calc(x,y,m):
    rmsd=0
    n=len(x)
    for i in range(n):
        rmsd+= (y[i]-x[i]*m)**2
    return math.sqrt(rmsd/n)


def rmsd_calc2(array,excol,eycol,m):
    rmsd=0
    n=len(array)
    for i in range(n):
        rmsd+= (float(array[i][eycol])-float(array[i][excol])*m)**2
    return math.sqrt(rmsd/n)


def pear(x,y):
    mean_x = mean(x)
    mean_y = mean(y)
    var_x = var(x)
    var_y = var(y)
    covar_xy = covar(x,y)
    #y=a*x+b, r is pearson coefficient
    out=[]
    out.append(covar_xy/var_x)                             #a
    out.append(-(covar_xy/var_x)*mean_x + mean_y)          #b
    out.append(covar_xy/math.sqrt(var_x)/math.sqrt(var_y)) #r
    return out

def pear2(array,col1,col2):
    x=[]
    y=[]
    for line in array:
        x.append(float(line[col1]))
        y.append(float(line[col2]))
    mean_x = mean(x)
    mean_y = mean(y)
    var_x = var(x)
    var_y = var(y)
    covar_xy = covar(x,y)
    #y=a*x+b, r is pearson coefficient
    out=[]
    out.append(covar_xy/var_x)                             #a
    out.append(-(covar_xy/var_x)*mean_x + mean_y)          #b
    out.append(covar_xy/math.sqrt(var_x)/math.sqrt(var_y)) #r
    return out


############################################################
