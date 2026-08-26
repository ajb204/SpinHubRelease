#!/usr/bin/python

#################################################################################
# Script for processing NMR data
# Uses NMR glue to read in a Fourier Transformed file in nmrPipe format
#
# A.Baldwin 21st December 2012
#################################################################################


import vpar
import numpy,sys,os,math
import nmrglue as ng
from scipy.optimize import leastsq
from baldwinStd import PathExists

#################################################################################

def AddIfNew(array,test):
    tick=0
    for i in range(len(array)):
        if(array[i]==test):
            tick=1
    if(tick==0):
        array.append(test)

def Find2Pow(gradl):
    for i in range(100):
        test=2**(i)
        if(test>len(gradl)):
            return test

def findnear_index(test,array):
    inear=0
    itest=math.fabs(float(array[0])-float(test))
    for i in range(len(array)):
        if(math.fabs(float(array[i])-float(test))<itest):
            itest=math.fabs(float(array[i])-float(test))
            inear=i
    return int(inear)



def FitFuncDiff(x,b):
    Deff=x[0]
    A0=x[1]
    grads=b[0]
    data=b[1]
    b=b[2]

    Ycalc=A0*numpy.exp(-grads*b*Deff)
    return Ycalc-data


def FitFuncDiff2(x,b):
    Deff=x[0]
    A0=x[1]
    Deff2=x[2]
    B0=x[3]

    grads=b[0]
    data=b[1]
    b=b[2]

    test='n'

    if(test=='y'):
        RetVal=[]
        for i in range(len(grads)):

            if(grads[i]<7800):
                RetVal.append((math.fabs(A0)*numpy.exp(-grads[i]*b*Deff)+math.fabs(B0)*numpy.exp(-grads[i]*b*Deff2)-data[i]))
        #print len(RetVal)
        return RetVal
    else:
        Ycalc=math.fabs(A0)*numpy.exp(-grads*b*Deff)+math.fabs(B0)*numpy.exp(-grads*b*Deff2)
        return Ycalc-data


def GetNoise(fac,data,index1,noiseMin,noiseMax):
    noise=fac*numpy.max(data[:,findnear_index(noiseMax,index1):findnear_index(noiseMin,index1)])
    print 'Noise:',noise
    return noise


def MakeProj(pth,data,index1,Size,noise):
    outy=open('out/proj.'+pth+'.out','w')
    for i in range(Size[1]):
        #outy.write('%f\t%f\t%f\n' % (index1[i],numpy.sum(data[:,i],axis=0),noise))
        outy.write('%f\t%f\t%f\n' % (index1[i],numpy.sum(data[0,i],axis=0),noise))
    outy.close()

def GetExptDetails(pth,Size):
    type=vpar.GetSpectrometerType(path=pth+'/')

    if(type=='var'):
        
        seqfil=vpar.GetParVarian(pth+'/procpar','n',('','seqfil'))[0].split('\"')[1]
        
        if(seqfil=="CT_N_hsqc_LED_lek_600_v2"):
            tau=float(vpar.GetParVarian(pth+'/procpar','n',('','Big_delta'))[0])
            gradt=float(vpar.GetParVarian(pth+'/procpar','n',('','gt2'))[0])
            gradl=vpar.GetParVarian(pth+'/procpar','n',('','gzlvl2'))
        elif(seqfil=="water_sLED_fm_v2_600"):
            tau=float(vpar.GetParVarian(pth+'/procpar','n',('','BigT'))[0])
            gradt=float(vpar.GetParVarian(pth+'/procpar','n',('','gt1'))[0])
            gradl=vpar.GetParVarian(pth+'/procpar','n',('','gzlvl1'))

        elif(seqfil=="hmqc_c13_600_methyl_diffusion_lek"):
            tau=float(vpar.GetParVarian(pth+'/procpar','n',('','bigT'))[0])
            gradt=2*float(vpar.GetParVarian(pth+'/procpar','n',('','gt5'))[0])
            gradl=vpar.GetParVarian(pth+'/procpar','n',('','gzlvl5'))           
        else:
            print 'cannot find pulse sequence',seqfil
            sys.exit(100)



    elif(type=='omega'):
        gradt=float(eval(vpar.GetOmegaVal(pth+'/test.par','gradt')[0].replace('ms','*1000')))
        tau=float(eval(vpar.GetOmegaVal(pth+'/test.par','tau')[0].replace('ms','*1000')))
        gradl=vpar.GetOmegaVal(pth+'/test.par','gradl')



        TIME_T2=float(vpar.GetParOmega(parfile,'n',('','time_T2'))[0])*1E-6


    print 'Gradient time:',gradt
    print 'BigT:',tau
    print 'Total grad pts:',len(gradl)

    try:
        grads=[]
        for i in range(Size[0]):
            grads.append(float(gradl[i%len(gradl)])**2.)

    except:
        gradMin=float(GetOmegaVal(pth+'/test.par','gradMin')[0])
        gradMax=float(GetOmegaVal(pth+'/test.par','gradMax')[0])
        nb=int(GetOmegaVal(pth+'/test.par','nb')[0])
        print 'GradMax: ',gradMax
        print 'GradMin: ',gradMin
        print 'nb:      ',nb
        grads=[]
        for i in range(nb):
            gradSq=gradMin*gradMin+i/(nb-1.)*(gradMax*gradMax-gradMin*gradMin)
            grads.append(gradSq)


    #adjust grads
    if(len(grads)!=Size[0]):
        gradNew=[]
        for i in range((Size[0])):
            gradNew.append(grads[i%(len(grads))])
    else:
            gradNew=grads
    

    gradNew=numpy.array(gradNew)


    return gradt,tau,gradl,gradNew


def CleanData(data,index1,ppmMax,ppmMin,Size):
    ppmMaxPt=findnear_index(ppmMax,index1)
    ppmMinPt=findnear_index(ppmMin,index1)
    ext=[]
    for j in range((Size[0])):
        test=numpy.sum(data[j,]*data[j,])
        if(test>0): #if sum of spectrum is not zero, keep it...
            ext.append(j)
            #print 'including'
        #else:
            #print 'skipping'
    strip= data[ext,ppmMaxPt:ppmMinPt]
    Size= strip.shape
    return Size,strip,ext


def FitPPM(pth,Size,strip,grads,bfac,ppmMax,noise,index1):
    print 'Analysing each ppm value...'
    outy=open('out/fitty.'+pth+'.out','w')
    outy2=open('out/fitty2.'+pth+'.out','w')
    ppmMaxPt=findnear_index(ppmMax,index1)

    for i in range(Size[1]):

        slice=numpy.array(strip[:,i])
        above= numpy.sum(slice>noise)
        if(above>3):
            x0=leastsq(FitFuncDiff,[1E-10,slice[0]],args=[grads,slice,bfac])
            outy.write('%f\t%e\t%f\n' % (index1[ppmMaxPt+i],x0[0][0],x0[0][1]))
            if(above>5): #if there are at least 5 points
                
                
                x0=leastsq(FitFuncDiff2,[x0[0][0],slice[0],x0[0][0]/5.,slice[0]/5.],args=[grads,slice,bfac])
                if(x0[0][0]>x0[0][3]):
                    outy2.write('%f\t%e\t%f\t%e\t%f\n' % (index1[ppmMaxPt+i],x0[0][0],x0[0][1],x0[0][2],x0[0][3]))
                else:
                    outy2.write('%f\t%e\t%f\t%e\t%f\n' % (index1[ppmMaxPt+i],x0[0][2],x0[0][3],x0[0][0],x0[0][1]))
    outy.close()





def FitInt(IntMin,IntMax,bfac,data,index1,ext,grads,pth):
    print 'Integrating between ',IntMin,'and',IntMax

    strip= data[ext,findnear_index(IntMax,index1):findnear_index(IntMin,index1)]
    integ= numpy.sum(strip,axis=1) #sum over range

    x0=leastsq(FitFuncDiff,[1E-10,integ[0]],args=[grads,integ,bfac])


    outy=open('out/sum.'+pth+'.norm.out','w')
    outy.write('# Diffusion: %e\n' % (x0[0][0]))
    outy.write('# A0: %e\n' % (x0[0][1]))
    for i in range(len(integ)):
        Xval=bfac*grads[i]
        outy.write('%s\t%f\n' % (Xval,integ[i]/x0[0][1]))
    outy.write('\n\n')
    gradmin=numpy.min(grads)
    gradmax=numpy.max(grads)
    for i in range((100)):
        Xval=bfac*(gradmin+i/(100.-1)*(gradmax-gradmin))
        outy.write('%s\t%f\n' % (Xval,math.exp(-Xval*x0[0][0])))        
    outy.close()

    outy=open('out/sum.'+pth+'.out','w')
    outy.write('# Diffusion: %e\n' % (x0[0][0]))
    outy.write('# A0: %e\n' % (x0[0][1]))
    for i in range(len(integ)):
        Xval=bfac*grads[i]
        outy.write('%s\t%f\n' % (Xval,integ[i]))
    outy.write('\n\n')
    gradmin=numpy.min(grads)
    gradmax=numpy.max(grads)
    for i in range((100)):
        Xval=bfac*(gradmin+i/(100.-1)*(gradmax-gradmin))
        outy.write('%s\t%f\n' % (Xval,x0[0][1]*math.exp(-Xval*x0[0][0])))        


    x1=leastsq(FitFuncDiff2,[x0[0][0],integ[0],x0[0][0]/5.,integ[0]/5.],args=[grads,integ,bfac])

    outy=open('out/sum2.'+pth+'.norm.out','w')
    outy.write('# Diffusion: %e\n' % (x1[0][0]))
    outy.write('# A0: %e\n' % (x1[0][1]))
    outy.write('# Diffusion2: %e\n' % (x1[0][2]))
    outy.write('# A0_2: %e\n' % (x1[0][3]))
    for i in range(len(integ)):
        Xval=bfac*grads[i]
        outy.write('%s\t%f\n' % (Xval,integ[i]/x0[0][1]))
    outy.write('\n\n')
    gradmin=numpy.min(grads)
    gradmax=numpy.max(grads)
    for i in range((100)):
        Xval=bfac*(gradmin+i/(100.-1)*(gradmax-gradmin))
        outy.write('%s\t%f\n' % (Xval,(x1[0][1]*math.exp(-Xval*x1[0][0])+x1[0][3]*math.exp(-Xval*x1[0][2]))/(x1[0][1]+x1[0][3]) ))        
    outy.close()

    outy=open('out/sum2.'+pth+'.out','w')
    outy.write('# Diffusion: %e\n' % (x1[0][0]))
    outy.write('# A0: %e\n' % (x1[0][1]))
    outy.write('# Diffusion2: %e\n' % (x1[0][2]))
    outy.write('# A0_2: %e\n' % (x1[0][3]))
    for i in range(len(integ)):
        Xval=bfac*grads[i]
        outy.write('%s\t%f\n' % (Xval,integ[i]))
    outy.write('\n\n')
    gradmin=numpy.min(grads)
    gradmax=numpy.max(grads)
    for i in range((100)):
        Xval=bfac*(gradmin+i/(100.-1)*(gradmax-gradmin))
        outy.write('%s\t%f\n' % (Xval,x1[0][1]*math.exp(-Xval*x1[0][0])+x1[0][3]*math.exp(-Xval*x1[0][2])


))        

    outy.close()


def ReadDiffData(infile):
    print 'Importing file ',infile
    dic,data=ng.pipe.read(infile)
    Size=data.shape
    uc0 = ng.pipe.make_uc(dic,data,dim=0)
    uc1 = ng.pipe.make_uc(dic,data,dim=1)
    print "Spectrum dimensions (pts): ",Size   #print the spectral dimensions
    print "dimension 0 limits (ppm): ", uc0.ppm(0), uc0.ppm(Size[0]-1)  #carbon (max/min)
    print "dimension 1 limits (ppm): ", uc1.ppm(0), uc1.ppm(Size[1]-1)  #direct proton (max/min)

    index0=[]#make index of carbon chemical shifts for index 0
    for i in range((Size[0])):
        index0.append((uc0.ppm(0)-i*(-uc0.ppm(Size[0]-1)+uc0.ppm(0))/(Size[0]-1)))
    index1=[]#make index of carbon chemical shifts for index 1
    for i in range((Size[1])):
        index1.append((uc1.ppm(0)-i*(-uc1.ppm(Size[1]-1)+uc1.ppm(0))/(Size[1]-1)))
    return data,Size,index0,index1


    


def DiffData(pth,noiseRng,noiseFac,extractRng,intRng):
    data,Size,index0,index1=ReadDiffData(pth+'/test.ft2')


    noise=GetNoise(noiseFac,data,index1,noiseRng[0],noiseRng[1])
    MakeProj(pth,data,index1,Size,noise)
    Size,strip,ext=CleanData(data,index1,extractRng[0],extractRng[1],Size)

    gradt,tau,gradl,grads=GetExptDetails(pth,Size)

    type=vpar.GetSpectrometerType(path=pth+'/')
    if(type=='var'):
        #normalise constants
#        yg=267.513  # *10**6 rad s-1 T-1
        yg=2.67513E4  # *rad s-1 G-1
        gmax= 0.002   #DAC to Gcm-1
        bfac=(yg*gmax*(gradt))**2*(tau-(gradt)/3)

    elif(type=='omega'):
        #normalise constants
        yg=267.513  # *10**6 rad s-1 T-1
        gmax= 60.   #G cm-1
        bfac=(yg*gmax*(gradt*1E-6))**2*(tau*1E-6-(gradt*1E-6)/3)

    FitPPM(pth,Size,strip,grads,bfac,extractRng[0],noise,index1) #fit each ppm
    FitInt(intRng[0],intRng[1],bfac,data,index1,ext,grads,pth)  #integrate over range
    return noise


def MakeGnuplot(pth,extractRng,integrate,specRng,noise):

    gnu=open('figs/gnu.gp','w')
    gnu.write('set term post eps enh color solid\n')
    gnu.write('set output \'figs/norm.eps\'\n')
    gnu.write('set title \'Normalised decay curve\'\n')
    gnu.write('set size square\n')
    gnu.write('set ylabel \'I/I_0\'\n')
    gnu.write('set xlabel \' ({/Symbol gd}G)^2({/Symbol D-d/3})/10^6\'\n')
    gnu.write('plot \\\n')
    for i in range(len(pth)):
        if(i!=0):
            gnu.write(',')
        gnu.write('\'out/sum.%s.norm.out\' i 0 u ($1/1E6):2 noti lt %i,\'\' i 1 u ($1/1E6):2 ti \'%s\' w li lt %i' % (pth[i],i+1,pth[i]+' '+str(integrate[0])+' to '+str(integrate[1])+' ppm',i+1))


    gnu.write('\n')
    gnu.write('reset\n')
    gnu.write('set term post eps enh color solid\n')
    gnu.write('set output \'figs/lin.eps\'\n')
    gnu.write('set title \'Normalised decay curve\'\n')
    gnu.write('set size square\n')
    gnu.write('set logscale y\n')
    gnu.write('set ylabel \'I/I_0\'\n')
    gnu.write('set xlabel \' ({/Symbol gd}G)^2({/Symbol D-d/3})/10^6\'\n')
    gnu.write('plot \\\n')
    for i in range(len(pth)):
        if(i!=0):
            gnu.write(',')
        gnu.write('\'out/sum.%s.out\' i 0 u ($1/1E6):2 noti lt %i,\'\' i 1 u ($1/1E6):2 ti \'%s\' w li lt %i' % (pth[i],i+1,pth[i]+' '+str(integrate[0])+' to '+str(integrate[1])+' ppm',i+1))

    for i in range(len(pth)):
        gnu.write(',')
        gnu.write('\'out/sum2.%s.out\' i 0 u ($1/1E6):2 noti lt %i,\'\' i 1 u ($1/1E6):2 ti \'%s\' w li lt %i' % (pth[i],i+1,pth[i]+' '+str(integrate[0])+' to '+str(integrate[1])+' ppm',i+1))


    for i in range(len(noise)):
        gnu.write(',\'\' i 0 u ($1/1E6):(%e) noti' % (noise[i]))

    gnu.write('\n')
    gnu.write('reset\n')
    gnu.write('set term post eps enh color solid\n')
    gnu.write('set output \'figs/diff.eps\'\n')
    gnu.write('set title \'Diffusion coefficients\'\n')
    gnu.write('set size square\n')
    gnu.write('set logscale y\n')
    gnu.write('set ylabel \'D_{eff} (cm^2s^{-1})\'\n')
    gnu.write('set xlabel \' {/Symbol d} (ppm)\'\n')
    gnu.write('set xrange[%s:%s]\n' % (str(extractRng[0]),str(extractRng[1])))
    gnu.write('set format y "10^{%T}"\n')
    gnu.write('plot \\\n')
    for i in range(len(pth)):
        if(i!=0):
            gnu.write(',')
        gnu.write('\'out/fitty.%s.out\' i 0 u 1:2 ti \'%s\' lt %i' % (pth[i],pth[i],i+1))




    gnu.write('\n')
    gnu.write('reset\n')
    gnu.write('set term post eps enh color solid\n')
    gnu.write('set title \'Projected spectrum\'\n')
    gnu.write('set size square\n')
    gnu.write('set key outside\n')
    gnu.write('set ylabel \'I\'\n')
    gnu.write('set xlabel \' {/Symbol d} (ppm)\'\n')
    gnu.write('set xrange[%s:%s]\n' % (str(extractRng[0]),str(extractRng[1])))
    for i in range(len(pth)):
        gnu.write('set output \'figs/proj.%s.eps\'\n'% pth[i])
        gnu.write('set yrange [%s:%s]\n' % (str(specRng[i][0]),str(specRng[i][1])))
        gnu.write('plot \'out/proj.%s.out\' u 1:2 ti \'projection\' w li lt %i,\'\' u 1:3 ti \'noise\' w li lt %i,\'out/fitty.%s.out\' u 1:3 ti \'%s\' w points  ps 0.3 pt 5 lt %i\n' % (pth[i],i+1,i+2,pth[i],'fitpts',i+3))

    gnu.close()

    os.system('gnuplot figs/gnu.gp')
    os.system('arraygraph.py 2 4 0 0 0 0 `ls figs/*.eps`')
    os.system('mv summary.pdf diffAnal.pdf')
    return

def AnalDiffDat(pth,noise,noiseFac,extract,integrate,specRng):



    PathExists(('out','figs'))
    

    noiselist=[]
    for i in range(len(pth)):
        noiselist.append(DiffData(pth[i],noise,noiseFac,extract,integrate))
    MakeGnuplot(pth,extract,integrate,specRng,noiselist)
    return noise
