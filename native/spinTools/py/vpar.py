#!/usr/bin/python

###########################################################
# Another version of vpar
# A.Baldwin 6th August 2012
#
# Need to check to make sure it can do long lists of things

import sys,os,string,numpy,math

def Find2Pow(gradl):
    for i in range(100):
        test=2**(i)
        if(test>=gradl):
            return test



def readfile(infile):
    peak=[]
    peakfile=open(infile,'r')
    for line in peakfile.readlines():
        linetosave=string.split(line)
        peak.append(linetosave)
    peakfile.close()
    return peak

#return water chemical shift in range 0-100oC
def WaterPPM(T):
    return 5.060 - 0.0122*T + (2.11E-5)*T**2.

# Based on Patrik Lundstrom 011126
#take water and sfrq, calc ppms of C and N
def Cshift(sfrq,dfrq,dfrq2,h1ppm,nuc='n'):
    C13_CONV=0.251449530
    N15_CONV=0.101329118
    P31_CONV=0.4048064954
    F19_CONV=0.9412866605363297

    if(nuc=='n'):


        sfrq0  = sfrq / (1.0 + h1ppm*1e-6);

        dfrq0 = sfrq0*C13_CONV;

        #dfrq0 = sfrq0*P31_CONV;
        dfrq20 = sfrq0*C13_CONV;

        c13ppm = (dfrq-dfrq0)/dfrq0*1e6;
        n15ppm = (dfrq2-dfrq20)/dfrq20*1e6;
        return c13ppm,n15ppm
    else:
        if(nuc=='P31'):
            conv=P31_CONV
        elif(nuc=='F19'):
            conv=F19_CONV
            print 'hello! I am F19',conv

        elif(nuc=='C13'):
            conv=C13_CONV
        elif(nuc=='N15'):
            conv=N15_CONV
        else:
            print 'could not find ',nuc
            sys.exit(100)
        print 'bastardo'
        sfrq0  = sfrq / (1.0 + h1ppm*1e-6);

        dfrq0 = sfrq0*conv;

        ppm = (dfrq-dfrq0)/dfrq0*1e6;
        return ppm
    



#take proton frequency and ppm,
#find carbon and nitrogen frequencies
#at specified ppms
def Cfrq(sfrq,dppm,dppm2,h1ppm):
    C13_CONV=0.251449530
    N15_CONV=0.101329118
    sfrq0  = sfrq / (1.0 + h1ppm*1e-6);

    dfrq0 = sfrq0*C13_CONV;
    dfrq20 = sfrq0*N15_CONV;

    dfrq=dfrq0*(dppm*1e-6+1.)
    dfrq2=dfrq20*(dppm2*1e-6+1.)

    return dfrq,dfrq2


#analyse either acqu and acqu2
def GetParBruk(infile,verb,argv):
    args=[]

    procpar=readfile(infile)
    for i in range(len(argv)-1):
        param=argv[i+1]
        tick=0
        for j in range(len(procpar)):
            test=procpar[j][0].split('##$')

            if(len(test)>1):
                test2=test[1].split('=')[0]
                if(test2==param):
                    if(verb=='y'):
                        sys.stdout.write('%s: %s\n' % (param,procpar[j][1]))
                    args.append(procpar[j][1])
                    tick=1
            else:
                #we have a line of zeros
                #is the previous line what we're after?
                test=procpar[j-1][0].split('##$')
                if(len(test)>1):
                    test2=test[1].split('=')[0]
                    for i in range(100):
                        parT=test2+str(i)
                        if(parT==param):
                            if(len(param.split(test2))>1):
                                if(verb=='y'):
                                    sys.stdout.write('Param %s found. Range: %s\n' % (param,procpar[j-1][1]))
                                #parameters are in rows in j,j+1,j+2...
                                go=0
                                cnt=0
                                while(go==0):
                                    if(i<len(procpar[j+cnt])):
                                        val=procpar[j+cnt][i]
                                        go=1
                                    else:
                                        i-=len(procpar[j+cnt])
                                        cnt+=1

                                args.append(val)
                                if(verb=='y'):
                                    sys.stdout.write('%s: %s\n' % (param,val))
                                tick=1
                #sys.exit(100)
    if(tick==0):
        if(verb=='y'):
            sys.stdout.write('Could not find param %s in %s\n' % (param,infile))
        return 'fail'
    else:
        return args




def GetParOmega(infile,verb,argv):
    args=[]
    procpar=readfile(infile)
    for i in range(len(argv)-1):
        param=argv[i+1]
        tick=0
        for j in range(len(procpar)):
            if(procpar[j][0]==param or procpar[j][1]==param):
                if(verb=='y'):
                    sys.stdout.write('%s: ' % (param))
                if(len(procpar[j])>2):
                    for k in range(len(procpar[j])-2):
                        if(procpar[j][k+1]!=param):
                            if(verb=='y'):
                                sys.stdout.write('%s ' % (procpar[j][k+1]))
                            args.append(procpar[j][k+1])
                else:
                    for k in range(len(procpar[j])-1):
                        if(verb=='y'):
                            sys.stdout.write('%s ' % (procpar[j][k+1]))
                        args.append(procpar[j][k+1])


                if(verb=='y'):
                    sys.stdout.write('\n')
                tick=1
    if(tick==0):
        if(verb=='y'):
            sys.stdout.write('Could not find param %s in %s\n' % (param,infile))
        return 'fail'
    else:
        return args






def GetParVarian(infile,verb,argv):
    args=[]
    procpar=readfile(infile)
    for i in range(len(argv)-1):
        param=argv[i+1]
        tick=0
        for j in range(len(procpar)):
            if(procpar[j][0]==param):
                if(verb=='y'):
                    sys.stdout.write('%s: %i argument' % (param,int(procpar[j+1][0])))
                if(int(procpar[j+1][0])>1):
                    sys.stdout.write('s')
                sys.stdout.write('\n')
                for k in range(int(procpar[j+1][0])):
                    try:
                        if(verb=='y'):
                            sys.stdout.write('%s ' % (procpar[j+1][k+1]))
                        args.append(procpar[j+1][k+1])
                    except:
                        if(verb=='y'):
                            sys.stdout.write('%s ' % (procpar[j+1+k][0]))
                        args.append(procpar[j+1+k][0])
                if(verb=='y'):
                    sys.stdout.write('\n')
                tick=1
    if(tick==0):
        if(verb=='y'):
            sys.stdout.write('Could not find param %s in procpar\n' % (param))
        return 'fail'
    else:
        return args


def GetParVarianFloat(infile,verb,argv):
    args=[]
    procpar=readfile(infile)
    for i in range(len(argv)-1):
        param=argv[i+1]
        tick=0
        for j in range(len(procpar)):
            if(procpar[j][0]==param):
                if(verb=='y'):
                    sys.stdout.write('%s: %i argument' % (param,int(procpar[j+1][0])))
                if(int(procpar[j+1][0])>1):
                    sys.stdout.write('s')
                sys.stdout.write('\n')
                for k in range(int(procpar[j+1][0])):
                    if(verb=='y'):
                        sys.stdout.write('%s ' % (procpar[j+1][k+1]))
                    args.append(float(procpar[j+1][k+1]))
                if(verb=='y'):
                    sys.stdout.write('\n')
                tick=1
    if(tick==0):
        if(verb=='y'):
            sys.stdout.write('Could not find param %s in procpar\n' % (param))
        return 'fail'
    else:
        return args



def GetOmegaVal(infile,param):
    if(GetParOmega(infile,'n',('',param))!='Fail'):
        return GetParOmega(infile,'n',('',param))


def GetBrukVal(infile,param):
    if(GetParBruk(infile,'n',('',param))!='Fail'):
        return GetParBruk(infile,'n',('',param))
    

    


def OmegaInfo(infile):
    sys.stdout.write('General information:\n\n')
    
    sys.stdout.write('Guessing number of dimensions...')

    na=int(GetOmegaVal(infile,'na')[0])
    nb=int(GetOmegaVal(infile,'nb')[0])
    nc=int(GetOmegaVal(infile,'nc')[0])

    np=int(GetOmegaVal(infile,'block_size')[0])
    if(nb==1 and nc==1):
        ndim=1
    elif(nb>1 and nc==1):
        ndi=2
    elif(nb>1 and nc>2):
        ndim=3
    sys.stdout.write('%i\n' % (ndim)) 
    sys.stdout.write('Direct dimension f1:\n')
    sfrq=float(GetOmegaVal(infile,'dim0_freq')[0])
    sw=float(GetOmegaVal(infile,'spec_width0')[0])
    sys.stdout.write('   sfrq:   %f MHz\n' % (sfrq))
    sys.stdout.write('   nt:     %i \n' % (na))
    sys.stdout.write('   at:     %f ms\n' % (float(GetOmegaVal(infile,'dwell_time0')[0])*1000.))
    sys.stdout.write('   np:     %i\n' % (np))
    sys.stdout.write('   sw:     %f Hz\n' % (float(GetOmegaVal(infile,'spec_width0')[0])))
    sys.stdout.write('   sw(P):  %f ppm\n' % (sw/sfrq))
    if(nb>1):
        sys.stdout.write('Inirect dimension f3:\n')
        dfrq=float(GetOmegaVal(infile,'f3_freq')[0])
        sw1=float(GetOmegaVal(infile,'sw1')[0])
        sys.stdout.write('   dfrq:   %f MHz\n' % (dfrq))
        sys.stdout.write('   ni:     %i \n' % (int(nb/2)))
        sys.stdout.write('   at:     %f ms\n' % (int(nb/2)*1/sw1*1000.))
        sys.stdout.write('   sw1:    %f Hz\n' % (sw1))
        sys.stdout.write('   sw1(P): %f ppm\n' % (sw1/dfrq))

    f1180_flg=int(GetOmegaVal(infile,'f1180_flg')[0])

    #NEED TO: get spectrometer frequency, nuclei and carrier
    ndim=2   #number of dimensions to convert
    nproc=2  #number of dimensions to fourier transform

    PipeParse('omega',infile.split('.par')[0]+'.bin',ndim,nproc,(2*np,nb),(np,nb/2),('Complex','Complex'),(sw,sw1),(sfrq,dfrq),('1H','13C'),f1180_flg)
    sys.stdout.write('\n')



def AddPipeLine(outy,lab,par,val,spa):
    outy.write(' %s%s ' % ('-'+lab+par.ljust(5),str(val).ljust(spa))) 

def EndPipeLine(outy):
    outy.write('\\\n')

def AddPipe(outy,ndim,axis,par,vals,spa):
    for i in range(ndim):
        AddPipeLine(outy,axis[i],par,vals[i],spa)
    EndPipeLine(outy)

def func(x,b):
    phi=x[0]
    data=b[0]
    shift=b[1]

    FT=numpy.fft.fft(data*complex(numpy.cos(phi),numpy.sin(phi)))
    
    #twin restraints:make the sum of the imag zero, and maximise the absorp
    b1=FT
    c=FT*b1*shift
    d=numpy.sum(c.imag)

    b1a=FT
    c1a=FT*b1a*shift
    d1a=numpy.sum(c1a.real)

    return 1.0*((d+d1a)*1E-6)



def AutoPhase(infid,ppmMin,ppmMax):
    import nmrglue as ng
    from scipy.optimize import leastsq
    dic,data = ng.pipe.read(infid)
    Size=data.shape
    
    uc0 = ng.pipe.make_uc(dic,data,dim=1)
    index0=[]#make index of carbon chemical shifts for index 0
    for i in range((Size[1])):
        index0.append((uc0.ppm(Size[1]-1)-i*(uc0.ppm(Size[1]-1)-uc0.ppm(0))/(Size[1]-1)))
    index0=numpy.array(index0)

    if(ppmMin=='*'):
        b3=index0>index0[0]
    else:
        b3=index0>ppmMin
    if(ppmMax=='*'):
        b2=index0<index0[len(index0)-1]
    else:
        b2=index0<ppmMax

    shift=b2*b3
    
    test=data[0]


    #outy=open('pants.out','w')
    for i in range(100):
        phi=i/100.*2*numpy.pi
        x0=leastsq(func,[phi,],args=[test,shift,])
        #x0=phi,
        phi=x0[0]
        val=func((x0[0],),(test,shift,))    
        #outy.write('%f\t%f\t%f\n' % (phi/numpy.pi*180.,x0[0],val))
        if(i==0):
            valTest=val
            phiVal=x0[0]
        if(val<valTest):
            valTest=val
            phiVal=x0[0]
    #outy.close()


    print 'Optimised phase angle:',phiVal/numpy.pi*180
    phi=phiVal


    print 'val:',func((phi,),(test,shift,))

    #phi=14.4/180.*numpy.pi
    #print 'val:',func((phi,),(test,shift,))


    FT=numpy.fft.fft(test*complex(numpy.cos(phi),numpy.sin(phi)))
    FT=numpy.fft.fftshift(FT)
    outy=open('test.out','w')
    for i in range(len(test)):
        outy.write('%f\t%f\n' % (index0[i],FT[i].real))
    outy.close()

    return phi/numpy.pi*180.


def PipeParse(type,infile,ndim,nproc,npT,np,mode,sw,frq,car,lab,f1180_flg,loop=0,ppmMin='*',ppmMax='*',collate_flg='n',f1neg='n'):    

    spa=10               #spacing in output file
    axis='x','y','z','a' #names of axes

    outy=open('./fid.test.com','w')
    outy.write('#!/bin/csh\n')
    if(loop==0):

        if(type=='omega'):
            outy.write('bin2pipe -in %s -ge -neg \\\n' % (infile))
        elif(type=='bruk'):
            outy.write('bruk2pipe -in %s  \\\n' % (infile))
            GRPDLY=float(GetBrukVal('acqus','GRPDLY')[0])
            DECIM=int(GetBrukVal('acqus','DECIM')[0])

            outy.write('-bad 0.0 -aswap -DMX -decim %i -dspfvs 20 -grpdly %f \\\n' % (DECIM,GRPDLY))
        elif(type=='var'):
            outy.write('var2pipe -in %s \\\n' % (infile))

    else: #loop over the four options if doing a spin-state selective expt
        outy.write('foreach f (0 1 2 3 )\n')
        if(type=='omega'):
            outy.write('bin2pipe -in %s.${f} -ge -neg \\\n' % (infile))
        elif(type=='bruk'):
            outy.write('bruk2pipe -in %s.${f}  \\\n' % (infile))
        elif(type=='var'):
            outy.write('var2pipe -in %s.${f} \\\n' % (infile))

    AddPipe(outy,ndim,axis,'N',npT,spa)
    AddPipe(outy,ndim,axis,'T',np,spa)
    AddPipe(outy,ndim,axis,'MODE',mode,spa)
    AddPipe(outy,ndim,axis,'SW',sw,spa)
    AddPipe(outy,ndim,axis,'OBS',frq,spa)
    AddPipe(outy,ndim,axis,'CAR',car,spa)
    AddPipe(outy,ndim,axis,'LAB',lab,spa)

    outy.write(' -ndim  %s -aq2D  %s \\\n' % (str(ndim).ljust(spa),'States'.ljust(spa)))
    if(collate_flg=='n'):
        if(loop==0):
            outy.write(' -out   test.fid -verb -ov\n')
        else:
            outy.write(' -out   test.fid.${f} -verb -ov\n')
            outy.write('end\n')
    else:
        outy.write('| nmrPipe -ov -verb -out test.fid\n') #spit into a giant fid
        

    #perform automatic phasing
    if(ppmMin=='n' or ppmMax=='n'):
        pha=0.
    else:
        try:
            pha=AutoPhase('test.fid',ppmMin=ppmMin,ppmMax=ppmMax)
        except:
            pha=0.

    outy=open('nmrproc.test.com','w')
    outy.write('#!/bin/csh\n')
    outy.write('nmrPipe -in test.fid \\\n')
    outy.write('#| nmrPipe -fn SOL                                 \\\n')
    if(type=='omega' or type=='var'): #for the omega, we multiply by 1 for the first point
        outy.write('#| nmrPipe  -fn EM  -lb 0.0 -c 1                    \\\n')
        outy.write('#| nmrPipe -fn SP  -off 0.5 -end 0.98 -pow 1 -c 1.0    \\\n')
        outy.write('| nmrPipe -fn GM -g1 2.0 -g2 10 -g3 0. -c 1.0           \\\n')
    else: #otherwise we multiply by 0.5
        outy.write('#| nmrPipe  -fn EM  -lb 0.0 -c 0.5                    \\\n')
        outy.write('#| nmrPipe -fn SP  -off 0.5 -end 0.98 -pow 1 -c 0.5    \\\n')
        outy.write('| nmrPipe -fn GM -g1 2.0 -g2 10 -g3 0. -c 0.5           \\\n')
    outy.write('| nmrPipe  -fn ZF -auto                            \\\n')
    outy.write('| nmrPipe  -fn FT -auto                            \\\n')
    outy.write('| nmrPipe  -fn PS -p0 %f -p1 0.00 -di -verb       \\\n' % pha)
    outy.write('#| nmrPipe -fn EXT -xn 3ppm -x1 -0.4ppm -sw         \\\n')

    if(nproc>1):
        outy.write('| nmrPipe  -fn TP                                  \\\n')
        outy.write('| nmrPipe  -fn LP -fb                              \\\n')
        if(f1180_flg==0):
            outy.write('#| nmrPipe  -fn EM  -lb 0.0 -c 0.5                   \\\n')
            outy.write('#| nmrPipe -fn SP  -off 0.5 -end 0.98 -pow 1 -c 0.5 \\\n')
            outy.write('| nmrPipe -fn GM -g1 2.0 -g2 10 -g3 0. -c 0.5       \\\n')
        else:
            outy.write('#| nmrPipe  -fn EM  -lb 0.0 -c 1.0                   \\\n')
            outy.write('#| nmrPipe -fn SP  -off 0.5 -end 0.98 -pow 1 -c 1.0  \\\n')
            outy.write('| nmrPipe -fn GM -g1 2.0 -g2 10 -g3 0. -c 1.0       \\\n')
        outy.write('| nmrPipe  -fn ZF -auto                            \\\n')
        if(f1neg=='n'):
            outy.write('| nmrPipe  -fn FT -auto                            \\\n')
        else:
            outy.write('| nmrPipe  -fn FT -neg                            \\\n')

        if(f1180_flg==0):
            outy.write('| nmrPipe  -fn PS -p0 0 -p1 0.00 -di -verb         \\\n')
        else:
            outy.write('| nmrPipe  -fn PS -p0 -90 -p1 180.00 -di -verb         \\\n')
        outy.write('#| nmrPipe -fn TP                                  \\\n')
    outy.write('#| nmrPipe -fn POLY -auto                          \\\n')
    outy.write('   -ov -out test.ft2\n')
    outy.close()
    return



def BrukInfo(infile1,infile2):
    sys.stdout.write('General information:\n\n')
    sys.stdout.write('Guessing number of dimensions...')

    nt=int(GetBrukVal(infile1,'NS')[0])
    np=int(GetBrukVal(infile1,'TD')[0])
    ni=int(GetBrukVal(infile2,'TD')[0])
    #nc=int(GetOmegaVal(infile,'nc')[0])
    #if(nb==1 and nc==1):
    #    ndim=1
    #elif(nb>1 and nc==1):
    #    ndim=2
    #elif(nb>1 and nc>2):
    #    ndim=3
    #sys.stdout.write('%i\n' % (ndim)) 
    sys.stdout.write('Direct dimension f1:\n')
    sfrq=float(GetBrukVal(infile1,'BF1')[0])
    sw=float(GetBrukVal(infile1,'SW_h')[0])
    sys.stdout.write('   sfrq:   %f MHz\n' % (sfrq))
    sys.stdout.write('   nt:     %i \n' % (nt))
    sys.stdout.write('   at:     %f ms\n' % (np/sw*1000.))
    sys.stdout.write('   np:     %i\n' % (int(GetBrukVal(infile1,'TD')[0])))
    sys.stdout.write('   sw:     %f Hz\n' % (sw))
    sys.stdout.write('   sw(P):  %f ppm\n' % (sw/sfrq))
    if(ni>1):
        sys.stdout.write('Inirect dimension f3:\n')
        dfrq=float(GetBrukVal(infile2,'BF1')[0])
        sw1=float(GetBrukVal(infile2,'SW_h')[0])
        sys.stdout.write('   dfrq:   %f MHz\n' % (dfrq))
        sys.stdout.write('   ni:     %i \n' % (ni))
        sys.stdout.write('   at:     %f ms\n' % (ni/sw1*1000.))
        sys.stdout.write('   sw1:    %f Hz\n' % (sw1))
        sys.stdout.write('   sw1(P): %f Hz\n' % (sw1/dfrq))
    sys.stdout.write('\n')
    #NEED TO: get spectrometer frequency, nuclei and carrier
    ndim=2   #number of dimensions to convert
    nproc=2  #number of dimensions to fourier transform
    
    zgopts=GetBrukVal(infile1,'ZGOPTNS')
    if(zgopts[0].split('F1180')>1):
        f1180_flg=1
    else:
        f1180_flg=0
    
    if(ndim==2):
        inny='ser'

    PipeParse('bruk',inny,ndim,nproc,(np,ni*2),(np/2,ni),('Complex','Complex'),(sw,sw1),(xcar,ycar),(sfrq,dfrq),('1H','13C'),f1180_flg)

    

def GetSpectrometerType(path='./'):
    if(os.path.exists(path+'acqu')==1 or os.path.exists('./acqu2')==1):
        sys.stdout.write('Found acqu: we are Bruker!\n')
        return 'bruk'
    elif(os.path.exists(path+'procpar')==1):
        sys.stdout.write('Found procpar: we are varian!\n')
        return 'var'

    else:
        sys.stdout.write('Neither bruker nor varian - guessing GE!\n')
        return 'omega'

def GetOmegaParFile():
    files=os.listdir('./')
    for i in range(len(files)):
        if(len(files[i].split('.par'))>1):
            sys.stdout.write('Found omega par file: %s\n' % files[i])
            return files[i]
    sys.stdout.write('Cannot find omega par file\n')
    sys.exit(100)

def GetParFile(path='./'):
    files=os.listdir(path)
    for i in range(len(files)):
        if(len(files[i].split('.par'))>1):
            return files[i]
    print 'Cannot find par file'
    return 0



def GetPar(path,pars):
    verb='y'

    parAdj=[]
    parAdj.append('')
    for i in range(len(pars)):
        parAdj.append(pars[i])

    if(os.path.exists(path+'procpar')==1):
        if(verb=='y'):
            sys.stdout.write('Found procpar: we are varian! Proceeding...\n')
        args=GetParVarian(path+'procpar',verb,parAdj)
    elif(os.path.exists(path+'acqu')==1 or os.path.exists('./acqu2')==1):
        if(verb=='y'):
            sys.stdout.write('Found acqu: we are Bruker! Proceeding...\n')
            sys.stdout.write('acqu:\n')
            BrukInfo(path+'acqu',path+'acqu2')
        args=GetParBruk('./acqu',verb,parAdj)
        if(verb=='y'):
            sys.stdout.write('acqu2:\n')
        args=GetParBruk(path+'acqu2',verb,parAdj)
    else:
        #we must be GE
        files=os.listdir(path)
        for i in range(len(files)):
            if(len(files[i].split('.par'))>1):

                if(verb=='y'):
                    sys.stdout.write('Found omega file: %s\n' % files[i])
                    OmegaInfo(files[i])
                args=GetParOmega(files[i],verb,parAdj)
    return args
    







if __name__ == '__main__':
    verb='y'
    if(len(sys.argv)==1):
        sys.stdout.write('USAGE: %s param1 param2 param3...\n' % (sys.argv[0]))
        sys.exit(100)

    if(os.path.exists('./procpar')==1):
        if(verb=='y'):
            sys.stdout.write('Found procpar: we are varian! Proceeding...\n')
        args=GetParVarian('./procpar',verb,sys.argv)
    elif(os.path.exists('./acqu')==1 or os.path.exists('./acqu2')==1):
        if(verb=='y'):
            sys.stdout.write('Found acqu: we are Bruker! Proceeding...\n')
            sys.stdout.write('acqu:\n')
            BrukInfo('./acqu','./acqu2')
        args=GetParBruk('./acqu',verb,sys.argv)
        if(verb=='y'):
            sys.stdout.write('acqu2:\n')
        args=GetParBruk('./acqu2',verb,sys.argv)
    else:
        #we must be GE
        files=os.listdir('./')
        for i in range(len(files)):
            if(len(files[i].split('.par'))>1):

                if(verb=='y'):
                    sys.stdout.write('Found omega file: %s\n' % files[i])
                    OmegaInfo(files[i])
                args=GetParOmega(files[i],verb,sys.argv)
        
    



