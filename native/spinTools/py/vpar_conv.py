#!/usr/bin/python

##############################################################
# Contains the conversion and processing functions for a 
# variety of standard experiments
#
# A.Baldwin 15th Feb 2014
#

import vpar,os,sys,reGlue,math


def ConvertSpec_HtoC(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','ncyc_cp')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(1)
    elif(array[1]=='phase'):
        mode=str(0)

    print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Complex','Real'),(sw,sw1,0.0),(sfrq,dfrq,0.0),(vpar.WaterPPM(temp),c13shift,0.0),('1H','13C','ncyc_cp'),f1180_flg)    

    os.system('csh fid.test.com')




def ConvertSpec_15Ncpmg(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','ncyc')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(1)
    elif(array[1]=='phase'):
        mode=str(0)

    print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

#    os.system('re-shuffle.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')

    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Rance-Kay','Real'),(sw,sw1,0.0),(sfrq,dfrq2,0.0),(vpar.WaterPPM(temp),n15shift,0.0),('1H','15N','ncyc'),f1180_flg,collate_flg='y',f1neg='y')    

    os.system('csh fid.test.com')




def Convertspec_15NHT2(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','dly_T2')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(1)
    elif(array[1]=='phase'):
        mode=str(0)

    print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

#    os.system('re-shuffle.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    

   # print 'RelaxFix.out ', str(np), ' '+str(ni), ' ', str(nz), ' ', mode, ' fid.final fid'

#    os.system('./RelaxFix.out '+str(np)+' fid fid.final')
 
    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Rance-Kay','Real'),(sw,sw1,0.0),(sfrq,dfrq2,0.0),(vpar.WaterPPM(temp),n15shift,0.0),('1H','15N','ncyc'),f1180_flg,collate_flg='y',f1neg='n')    

    os.system('csh fid.test.com')



def Convertspec_15NHNOE(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','ncyc')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(1)
    elif(array[1]=='phase'):
        mode=str(0)

    print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

#    os.system('re-shuffle.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    

   # print 'RelaxFix.out ', str(np), ' '+str(ni), ' ', str(nz), ' ', mode, ' fid.final fid'

#    os.system('./RelaxFix.out '+str(np)+' fid fid.final')
 
    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Rance-Kay','Real'),(sw,sw1,0.0),(sfrq,dfrq2,0.0),(vpar.WaterPPM(temp),n15shift,0.0),('1H','15N','ncyc'),f1180_flg,collate_flg='y',f1neg='n')    

    os.system('csh fid.test.com')


def Convertspec_15NHT2b(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','ncyc')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(1)
    elif(array[1]=='phase'):
        mode=str(0)

    print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

#    os.system('re-shuffle.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    

   # print 'RelaxFix.out ', str(np), ' '+str(ni), ' ', str(nz), ' ', mode, ' fid.final fid'

#    os.system('./RelaxFix.out '+str(np)+' fid fid.final')
 
    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Rance-Kay','Real'),(sw,sw1,0.0),(sfrq,dfrq2,0.0),(vpar.WaterPPM(temp),n15shift,0.0),('1H','15N','ncyc'),f1180_flg,collate_flg='y',f1neg='n')    

    os.system('csh fid.test.com')



def Convertspec_15NHT1(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','ncyc')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(1)
    elif(array[1]=='phase'):
        mode=str(0)

    print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

#    os.system('re-shuffle.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    

   # print 'RelaxFix.out ', str(np), ' '+str(ni), ' ', str(nz), ' ', mode, ' fid.final fid'

#    os.system('./RelaxFix.out '+str(np)+' fid fid.final')
 
    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Rance-Kay','Real'),(sw,sw1,0.0),(sfrq,dfrq2,0.0),(vpar.WaterPPM(temp),n15shift,0.0),('1H','15N','ncyc'),f1180_flg,collate_flg='y',f1neg='n')    

    os.system('csh fid.test.com')
	


def Convertspec_15NHT1unenhanced(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','ncyc')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(1)
    elif(array[1]=='phase'):
        mode=str(0)

    print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

#    os.system('re-shuffle.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    

   # print 'RelaxFix.out ', str(np), ' '+str(ni), ' ', str(nz), ' ', mode, ' fid.final fid'

#    os.system('./RelaxFix.out '+str(np)+' fid fid.final')
 
    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Complex','Real'),(sw,sw1,0.0),(sfrq,dfrq2,0.0),(vpar.WaterPPM(temp),n15shift,0.0),('1H','15N','ncyc'),f1180_flg,collate_flg='y',f1neg='n')    

    os.system('csh fid.test.com')
	





def ConvertSpec_15NcpmgTrosyAntiTrosy(ni=0,ppmMin='*',ppmMax='*'):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','ncyc')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(1)
    elif(array[1]=='phase'):
        mode=str(0)

    print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

#    os.system('re-shuffle.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')

    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Rance-Kay','Real'),(sw,sw1,0.0),(sfrq,dfrq2,0.0),(vpar.WaterPPM(temp),c13shift,0.0),('1H','15N','ncyc'),f1180_flg,collate_flg='y',f1neg='y')    

    os.system('csh fid.test.com')




def ConvertSpec_15Ncpmg_omeg(infile,ni=0):

    na=int(vpar.GetOmegaVal(infile,'na')[0])
    nb=int(vpar.GetOmegaVal(infile,'nb')[0])
    nc=int(vpar.GetOmegaVal(infile,'nc')[0])

    np=int(vpar.GetOmegaVal(infile,'block_size')[0])

    if(nc==1): #if we're only acquiring one plane
        ndim=2
        sys.stdout.write('%i\n' % (ndim)) 
        sys.stdout.write('Direct dimension f1:\n')
        sfrq=float(vpar.GetOmegaVal(infile,'dim0_freq')[0])
        sw=float(vpar.GetOmegaVal(infile,'spec_width0')[0])
        sys.stdout.write('   sfrq:   %f MHz\n' % (sfrq))
        sys.stdout.write('   nt:     %i \n' % (na))
        sys.stdout.write('   at:     %f ms\n' % (float(vpar.GetOmegaVal(infile,'dwell_time0')[0])*1000.))
        sys.stdout.write('   np:     %i\n' % (np))
        sys.stdout.write('   sw:     %f Hz\n' % (float(vpar.GetOmegaVal(infile,'spec_width0')[0])))
        sys.stdout.write('   sw(P):  %f ppm\n' % (sw/sfrq))
        
        sys.stdout.write('Inirect dimension f2:\n')
        dfrq=float(vpar.GetOmegaVal(infile,'f2_freq')[0])
        sw1=float(vpar.GetOmegaVal(infile,'sw1')[0])
        sys.stdout.write('   dfrq:   %f MHz\n' % (dfrq))
        sys.stdout.write('   ni:     %i \n' % (int(nb/2)))
        sys.stdout.write('   at:     %f ms\n' % (int(nb/2)*1/sw1*1000.))
        sys.stdout.write('   sw1:    %f Hz\n' % (sw1))
        sys.stdout.write('   sw1(P): %f ppm\n' % (sw1/dfrq))
        
        print vpar.GetOmegaVal(infile,'f1180_flg')[0]
        f1180_flg=int(float(vpar.GetOmegaVal(infile,'f1180_flg')[0])*1.0)
        
        mode=str(1)
        
    #os.system('RelaxFix.out '+str(np)+' '+str(nc/2)+' '+str(nb)+' '+mode+' fid.final')
        
    #NEED TO: get spectrometer frequency, nuclei and carrier
        nproc=2  #number of dimensions to fourier transform
        
        temp=float(vpar.GetOmegaVal(infile,'temp')[0])
        print 'water:',vpar.WaterPPM(temp)
        
        c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(temp))
        
        ndim=2 #first convert the whole thing as a pseudo 2D
        vpar.PipeParse('omega',infile.split('.par')[0]+'.bin',ndim,nproc,(2*np,nb),(np,nb*0.5),('Complex','Rance-Kay'),(sw,sw1),(sfrq,dfrq),(vpar.WaterPPM(temp),n15shift),('1H','15N'),f1180_flg) #extract the raw data as 2d
        
        print 'Running conversion script...'
        os.system('csh fid.test.com')
        
    #remove zero fids and reshuffle
    #manually change the header to convert to 3D


    else:
        ndim=3
        sys.stdout.write('%i\n' % (ndim)) 
        sys.stdout.write('Direct dimension f1:\n')
        sfrq=float(vpar.GetOmegaVal(infile,'dim0_freq')[0])
        sw=float(vpar.GetOmegaVal(infile,'spec_width0')[0])
        sys.stdout.write('   sfrq:   %f MHz\n' % (sfrq))
        sys.stdout.write('   nt:     %i \n' % (na))
        sys.stdout.write('   at:     %f ms\n' % (float(vpar.GetOmegaVal(infile,'dwell_time0')[0])*1000.))
        sys.stdout.write('   np:     %i\n' % (np))
        sys.stdout.write('   sw:     %f Hz\n' % (float(vpar.GetOmegaVal(infile,'spec_width0')[0])))
        sys.stdout.write('   sw(P):  %f ppm\n' % (sw/sfrq))
        
        sys.stdout.write('Inirect dimension f2:\n')
        dfrq=float(vpar.GetOmegaVal(infile,'f2_freq')[0])
        sw1=float(vpar.GetOmegaVal(infile,'sw1')[0])
        sys.stdout.write('   dfrq:   %f MHz\n' % (dfrq))
        sys.stdout.write('   ni:     %i \n' % (int(nc/2)))
        sys.stdout.write('   at:     %f ms\n' % (int(nc/2)*1/sw1*1000.))
        sys.stdout.write('   sw1:    %f Hz\n' % (sw1))
        sys.stdout.write('   sw1(P): %f ppm\n' % (sw1/dfrq))
        
        print vpar.GetOmegaVal(infile,'f1180_flg')[0]
        f1180_flg=int(float(vpar.GetOmegaVal(infile,'f1180_flg')[0])*1.0)
        
        mode=str(1)
        
    #os.system('RelaxFix.out '+str(np)+' '+str(nc/2)+' '+str(nb)+' '+mode+' fid.final')
        
    #NEED TO: get spectrometer frequency, nuclei and carrier
        nproc=2  #number of dimensions to fourier transform
        
        temp=float(vpar.GetOmegaVal(infile,'temp')[0])
        print 'water:',vpar.WaterPPM(temp)
        
        c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(temp))
        
        nbP2=vpar.Find2Pow(nb) #round up because of stupid omega format
        ncP2=vpar.Find2Pow(nc) #round up because of stupid omega format
        
        ndim=2 #first convert the whole thing as a pseudo 2D
        vpar.PipeParse('omega',infile.split('.par')[0]+'.bin',ndim,nproc,(2*np,ncP2*nbP2),(np,ncP2*nbP2),('Complex','Complex'),(sw,sw1),(sfrq,dfrq),(vpar.WaterPPM(temp),n15shift),('1H','15N'),f1180_flg) #extract the raw data as 2d
        
        print 'Running conversion script...'
        os.system('csh fid.test.com')
        
    #remove zero fids and reshuffle
    #manually change the header to convert to 3D
        reGlue.reShuff('test.fid',nb,nc,RanceKay=1) 



def ConvertSpec_15Ncpmg_bruk(infile,ni=0):

    na=int(vpar.GetBrukVal('acqus','NS')[0])
    nb=int(vpar.GetBrukVal('acqu2s','TD')[0])
    nc=int(vpar.GetBrukVal('acqu3s','TD')[0])

    np=int(vpar.GetBrukVal('acqus','TD')[0])

    if(nc==1): #if we're only acquiring one plane
        ndim=2
        sys.stdout.write('%i\n' % (ndim)) 
        sys.stdout.write('Direct dimension f1:\n')
        
        sfrq=float(vpar.GetBrukVal('acqus','SFO1')[0])
        sw=float(vpar.GetBrukVal('acqus','SW_h')[0])

        sys.stdout.write('   sfrq:   %f MHz\n' % (sfrq))
        sys.stdout.write('   nt:     %i \n' % (na))
        sys.stdout.write('   at:     %f ms\n' % (np/sw))
        sys.stdout.write('   np:     %i\n' % (np))
        sys.stdout.write('   sw:     %f Hz\n' % (sw))
        sys.stdout.write('   sw(P):  %f ppm\n' % (sw/sfrq))
        
        sys.stdout.write('Inirect dimension f2:\n')
        dfrq=float(vpar.GetBrukVal('acqu2','SFO1')[0])
        sw1=float(vpar.GetBrukVal('acqu2','SW_h')[0])

        sys.stdout.write('   dfrq:   %f MHz\n' % (dfrq))
        sys.stdout.write('   ni:     %i \n' % (int(nb/2)))
        sys.stdout.write('   at:     %f ms\n' % (int(nb/2)*1/sw1*1000.))
        sys.stdout.write('   sw1:    %f Hz\n' % (sw1))
        sys.stdout.write('   sw1(P): %f ppm\n' % (sw1/dfrq))
        
        #print vpar.GetOmegaVal(infile,'f1180_flg')[0]
        #f1180_flg=int(float(vpar.GetOmegaVal(infile,'f1180_flg')[0])*1.0)
        f1180_flg=1
        mode=str(1)
        
    #os.system('RelaxFix.out '+str(np)+' '+str(nc/2)+' '+str(nb)+' '+mode+' fid.final')
        
    #NEED TO: get spectrometer frequency, nuclei and carrier
        nproc=2  #number of dimensions to fourier transform
        
        temp=float(vpar.GetBrukVal('acqu','TEMP')[0])
        print 'water:',vpar.WaterPPM(temp)
        
        c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(temp))
        
        ndim=2 #first convert the whole thing as a pseudo 2D
        vpar.PipeParse('bruk','ser',ndim,nproc,(2*np,nb),(np,nb*0.5),('Complex','Rance-Kay'),(sw,sw1),(sfrq,dfrq),(vpar.WaterPPM(temp),n15shift),('1H','15N'),f1180_flg) #extract the raw data as 2d
        
        print 'Running conversion script...'
        os.system('csh fid.test.com')
        
    #remove zero fids and reshuffle
    #manually change the header to convert to 3D


    else:
        ndim=3
        sys.stdout.write('%i\n' % (ndim)) 
        sys.stdout.write('Direct dimension f1:\n')
        sfrq=float(vpar.GetBrukVal('acqus','SFO1')[0])
        sw=float(vpar.GetBrukVal('acqus','SW_h')[0])
        sys.stdout.write('   sfrq:   %f MHz\n' % (sfrq))
        sys.stdout.write('   nt:     %i \n' % (na))
        sys.stdout.write('   at:     %f ms\n' % (np/sw))
        sys.stdout.write('   np:     %i\n' % (np))
        sys.stdout.write('   sw:     %f Hz\n' % (sw))
        sys.stdout.write('   sw(P):  %f ppm\n' % (sw/sfrq))
        
        sys.stdout.write('Inirect dimension f2:\n')
        dfrq=float(vpar.GetBrukVal('acqu2s','SFO1')[0])
        sw1=float(vpar.GetBrukVal('acqu2s','SW_h')[0])

        sys.stdout.write('   dfrq:   %f MHz\n' % (dfrq))
        sys.stdout.write('   ni:     %i \n' % (int(nc/2)))
        sys.stdout.write('   at:     %f ms\n' % (int(nc/2)*1/sw1*1000.))
        sys.stdout.write('   sw1:    %f Hz\n' % (sw1))
        sys.stdout.write('   sw1(P): %f ppm\n' % (sw1/dfrq))
        
        #print vpar.GetOmegaVal(infile,'f1180_flg')[0]
        #f1180_flg=int(float(vpar.GetOmegaVal(infile,'f1180_flg')[0])*1.0)
        f1180_flg=1
        mode=str(1)
        
    #os.system('RelaxFix.out '+str(np)+' '+str(nc/2)+' '+str(nb)+' '+mode+' fid.final')
        
    #NEED TO: get spectrometer frequency, nuclei and carrier
        nproc=2  #number of dimensions to fourier transform

        try:
            temp=float(vpar.GetBrukVal('acqu','TEMP')[0])
        except:
            temp=30
        
        print 'water:',vpar.WaterPPM(temp)
        
        cnst18=float(vpar.GetBrukVal('acqus','CNST18')[0])
        cnst19=float(vpar.GetBrukVal('acqus','CNST19')[0])

        c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(temp))
        h1shift=vpar.WaterPPM(temp)-cnst18+cnst19
        
        npAdj=BrukFidAdjust(np)

        

        ndim=2 #first convert the whole thing as a pseudo 2D
#        vpar.PipeParse('bruk','ser',ndim,nproc,(npAdj,nc*nb),(np/2,nc*nb),('DQD','Complex'),(sw,sw1),(sfrq,dfrq),(vpar.WaterPPM(temp),n15shift),('1H','15N'),f1180_flg) #extract the raw data as 2d


        vpar.PipeParse('bruk','ser',ndim,nproc,(npAdj,nc*nb),(np/2,nc*nb),('DQD','Complex'),(sw,sw1),(sfrq,dfrq),(h1shift,n15shift),('1H','15N'),f1180_flg) #extract the raw data as 2d
        
        print 'Running conversion script...'
        os.system('csh fid.test.com')
        
    #remove zero fids and reshuffle
    #manually change the header to convert to 3D
        #reGlue.reShuff('test.fid',nb,nc,RanceKay=1) 

        reGlue.reShuff('test.fid',nc,nb,RanceKay=1) 

        









def ConvertSpec_15N_LED(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','gzlvl2')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    print array
    if(int(ni)==1): #if we're a 1D only expeirment...
        
        print 'ni:',ni,'np:',np,'nz:',nz

        sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
        sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])
        
        sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
        dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
        dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
        f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]
        
        temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
        print 'water:',vpar.WaterPPM(temp)
        
        c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))
        
        vpar.PipeParse('var','fid',2,1,(np,nz),(np/2,nz),('Complex','Real'),(sw,0.0),(sfrq,0.0),(vpar.WaterPPM(temp),0.0),('1H','gzlvl'),f1180_flg)    

    else:


        if(array[0]=='phase'):
            mode=str(1)
        elif(array[len(array)-1]=='phase'):
            mode=str(0)

        print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

        sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
        sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])
        
        sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
        dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
        dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
        f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]
        
        temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
        print 'water:',vpar.WaterPPM(temp)
        
        c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))
        
        ip_flg=vpar.GetParVarian('./procpar','n',('','IP_flg'))
        if(len(ip_flg)==2):
            os.system('re-shuffle2.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.fix fid')
        
            os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid.fix')
        else:
            os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')

        vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Complex','Real'),(sw,sw1,0.0),(sfrq,dfrq2,0.0),(vpar.WaterPPM(temp),n15shift,0.0),('1H','15N','gzlvl'),f1180_flg)    

        

    os.system('csh fid.test.com')




def ConvertSpec_13C_diff(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','gzlvl5')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    print array
    if(int(ni)==1): #if we're a 1D only expeirment...
        
        print 'ni:',ni,'np:',np,'nz:',nz

        sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
        sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])
        
        sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
        dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
        dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
        f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]
        
        temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
        print 'water:',vpar.WaterPPM(temp)
        
        c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))
        
        vpar.PipeParse('var','fid',2,1,(np,nz),(np/2,nz),('Complex','Real'),(sw,0.0),(sfrq,0.0),(vpar.WaterPPM(temp),0.0),('1H','gzlvl'),f1180_flg)    

    else:


        if(array[0]=='phase'):
            mode=str(1)
        elif(array[len(array)-1]=='phase'):
            mode=str(0)

        print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

        sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
        sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])
        
        sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
        dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
        dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
        f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]
        
        temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
        print 'water:',vpar.WaterPPM(temp)
        
        c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))
        
        #os.system('re-shuffle2.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.fix fid')
        
        os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
        vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Complex','Real'),(sw,sw1,0.0),(sfrq,dfrq,0.0),(vpar.WaterPPM(temp),c13shift,0.0),('1H','13C','gzlvl'),f1180_flg)    

        

    os.system('csh fid.test.com')



def ConvertSpec_waterSLED(ni=0,ppmMin='*',ppmMax='*'):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','gzlvl1')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    print array
    if(int(ni)==1 or int(ni)==0): #if we're a 1D only expeirment...
        
        print 'ni:',ni,'np:',np,'nz:',nz

        sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
        sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])
        
        sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
        dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
        dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
        f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]
        
        temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
        print 'water:',vpar.WaterPPM(temp)
        
        c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))
        
        vpar.PipeParse('var','fid',2,1,(np,nz),(np/2,nz),('Complex','Real'),(sw,0.0),(sfrq,0.0),(vpar.WaterPPM(temp),0.0),('1H','gzlvl'),f1180_flg)    

    else:


        if(array[0]=='phase'):
            mode=str(1)
        elif(array[len(array)-1]=='phase'):
            mode=str(0)

        print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

        sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
        sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])
        
        sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
        dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
        dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
        f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]
        
        temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
        print 'water:',vpar.WaterPPM(temp)
        
        c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))
        
        #os.system('re-shuffle2.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.fix fid')
        
        os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
        vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Complex','Real'),(sw,sw1,0.0),(sfrq,dfrq,0.0),(vpar.WaterPPM(temp),c13shift,0.0),('1H','13C','gzlvl'),f1180_flg)    

        

    os.system('csh fid.test.com')







def ConvertSpec_chmqc(ni=0,ppmMin='*',ppmMax='*'):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])


    print 'ni:',ni,'np:',np

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    tof=float(vpar.GetParVarian('./procpar','n',('','tof'))[0])
    tof_me=float(vpar.GetParVarian('./procpar','n',('','tof_me'))[0])
    direct=(tof_me-tof)/sfrq+vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

    vpar.PipeParse('var','fid',2,2,(np,ni*2),(np/2,ni),('Complex','Complex'),(sw,sw1),(sfrq,dfrq),(direct,c13shift),('1H','13C'),f1180_flg,ppmMin=ppmMin,ppmMax=ppmMax)    

    os.system('csh fid.test.com')



def ConvertSpec_ct_chsqc(ni=0,ppmMin='*',ppmMax='*'):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])


    print 'ni:',ni,'np:',np

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    tof=float(vpar.GetParVarian('./procpar','n',('','tof'))[0])
    #tof_me=float(vpar.GetParVarian('./procpar','n',('','tof_me'))[0])
    #direct=(tof_me-tof)/sfrq+vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

    vpar.PipeParse('var','fid',2,2,(np,ni*2),(np/2,ni),('Complex','Complex'),(sw,sw1),(sfrq,dfrq),(vpar.WaterPPM(temp),c13shift),('1H','13C'),f1180_flg,ppmMin=ppmMin,ppmMax=ppmMax)    

    os.system('csh fid.test.com')



def ConvertSpec_nhsqc(ni=0,ppmMin='*',ppmMax='*'):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])


    print 'ni:',ni,'np:',np

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    tof=float(vpar.GetParVarian('./procpar','n',('','tof'))[0])
    #tof_me=float(vpar.GetParVarian('./procpar','n',('','tof_me'))[0])
    #direct=(tof_me-tof)/sfrq+vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

    if(ni==1):
        vpar.PipeParse('var','fid',1,1,(np,),(np/2,),('Complex',),(sw,),(sfrq,),(vpar.WaterPPM(temp),),('1H',),f1180_flg,ppmMin=ppmMin,ppmMax=ppmMax)    
    else:
        vpar.PipeParse('var','fid',2,2,(np,ni*2),(np/2,ni),('Complex','Rance-Kay'),(sw,sw1),(sfrq,dfrq2),(vpar.WaterPPM(temp),n15shift),('1H','15N'),f1180_flg,ppmMin=ppmMin,ppmMax=ppmMax)    

    os.system('csh fid.test.com')



def ConvertSpec_nhsqc(ni=0,ppmMin='*',ppmMax='*'):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])


    print 'ni:',ni,'np:',np

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    tof=float(vpar.GetParVarian('./procpar','n',('','tof'))[0])
    #tof_me=float(vpar.GetParVarian('./procpar','n',('','tof_me'))[0])
    #direct=(tof_me-tof)/sfrq+vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

    if(ni==1):
        vpar.PipeParse('var','fid',1,1,(np,),(np/2,),('Complex',),(sw,),(sfrq,),(vpar.WaterPPM(temp),),('1H',),f1180_flg,ppmMin=ppmMin,ppmMax=ppmMax)    
    else:
        vpar.PipeParse('var','fid',2,2,(np,ni*2),(np/2,ni),('Complex','Rance-Kay'),(sw,sw1),(sfrq,dfrq2),(vpar.WaterPPM(temp),n15shift),('1H','15N'),f1180_flg,ppmMin=ppmMin,ppmMax=ppmMax)    

    os.system('csh fid.test.com')


def ConvertSpec_noesy(ni=0,ppmMin='n',ppmMax='n'):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    nz=len(vpar.GetParVarian('./procpar','n',('','mix')))
    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(1)
    elif(array[1]=='phase'):
        mode=str(0)

    print vpar.GetParVarian('./procpar','n',('','seqfil')),'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

    tof=float(vpar.GetParVarian('./procpar','n',('','tof'))[0])
    tof2=float(vpar.GetParVarian('./procpar','n',('','tof2'))[0])
    #direct=(tof_me-tof)/sfrq+vpar.WaterPPM(temp)
    direct=vpar.WaterPPM(temp)
    indirect=(tof-tof2)/sfrq+vpar.WaterPPM(temp)
    print direct,indirect
    print 'Calling RelaxFix.out:', str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid'
    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    print 'done'
    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Complex','Real'),(sw,sw1,0.0),(sfrq,sfrq,0.0),(direct,indirect,0.0),('1H1','1H2','mix'),f1180_flg,ppmMin=ppmMin,ppmMax=ppmMax)    
    print 'done pipeparse'
    os.system('csh fid.test.com')


def ConvertSpec_HtoC_omeg(infile,ni=0):

    na=int(vpar.GetOmegaVal(infile,'na')[0])
    nb=int(vpar.GetOmegaVal(infile,'nb')[0])
    nc=int(vpar.GetOmegaVal(infile,'nc')[0])

    np=int(vpar.GetOmegaVal(infile,'block_size')[0])


    ndim=3
    sys.stdout.write('%i\n' % (ndim)) 
    sys.stdout.write('Direct dimension f1:\n')
    sfrq=float(vpar.GetOmegaVal(infile,'dim0_freq')[0])
    sw=float(vpar.GetOmegaVal(infile,'spec_width0')[0])
    sys.stdout.write('   sfrq:   %f MHz\n' % (sfrq))
    sys.stdout.write('   nt:     %i \n' % (na))
    sys.stdout.write('   at:     %f ms\n' % (float(vpar.GetOmegaVal(infile,'dwell_time0')[0])*1000.))
    sys.stdout.write('   np:     %i\n' % (np))
    sys.stdout.write('   sw:     %f Hz\n' % (float(vpar.GetOmegaVal(infile,'spec_width0')[0])))
    sys.stdout.write('   sw(P):  %f ppm\n' % (sw/sfrq))

    sys.stdout.write('Inirect dimension f3:\n')
    dfrq=float(vpar.GetOmegaVal(infile,'f3_freq')[0])
    sw1=float(vpar.GetOmegaVal(infile,'sw1')[0])
    sys.stdout.write('   dfrq:   %f MHz\n' % (dfrq))
    sys.stdout.write('   ni:     %i \n' % (int(nb/2)))
    sys.stdout.write('   at:     %f ms\n' % (int(nb/2)*1/sw1*1000.))
    sys.stdout.write('   sw1:    %f Hz\n' % (sw1))
    sys.stdout.write('   sw1(P): %f ppm\n' % (sw1/dfrq))

    print vpar.GetOmegaVal(infile,'f1180_flg')[0]
    f1180_flg=int(float(vpar.GetOmegaVal(infile,'f1180_flg')[0])*1.0)


    mode=str(1)

    #os.system('RelaxFix.out '+str(np)+' '+str(nc/2)+' '+str(nb)+' '+mode+' fid.final')

    #NEED TO: get spectrometer frequency, nuclei and carrier
    nproc=2  #number of dimensions to fourier transform

    temp=float(vpar.GetOmegaVal(infile,'temp')[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(temp))

    nbP2=vpar.Find2Pow(nb) #round up because of stupid omega format
    ncP2=vpar.Find2Pow(nc) #round up because of stupid omega format

    ndim=2 #first convert the whole thing as a pseudo 2D
    vpar.PipeParse('omega',infile.split('.par')[0]+'.bin',ndim,nproc,(2*np,ncP2*nbP2),(np,ncP2*nbP2),('Complex','Complex'),(sw,sw1),(sfrq,dfrq),(vpar.WaterPPM(temp),c13shift),('1H','13C'),f1180_flg) #extract the raw data as 2d

    print 'Running conversion script...'
    os.system('csh fid.test.com')

    #remove zero fids and reshuffle
    #manually change the header to convert to 3D
    reGlue.reShuff('test.fid',nb,nc) 
    


def ConvertSpec_chmqc_omeg(infile,ni=0):

    na=int(vpar.GetOmegaVal(infile,'na')[0])
    nb=int(vpar.GetOmegaVal(infile,'nb')[0])
    nc=int(vpar.GetOmegaVal(infile,'nc')[0])

    np=int(vpar.GetOmegaVal(infile,'block_size')[0])


    ndim=2
    sys.stdout.write('%i\n' % (ndim)) 
    sys.stdout.write('Direct dimension f1:\n')
    sfrq=float(vpar.GetOmegaVal(infile,'dim0_freq')[0])
    sw=float(vpar.GetOmegaVal(infile,'spec_width0')[0])
    sys.stdout.write('   sfrq:   %f MHz\n' % (sfrq))
    sys.stdout.write('   nt:     %i \n' % (na))
    sys.stdout.write('   at:     %f ms\n' % (float(vpar.GetOmegaVal(infile,'dwell_time0')[0])*1000.))
    sys.stdout.write('   np:     %i\n' % (np))
    sys.stdout.write('   sw:     %f Hz\n' % (float(vpar.GetOmegaVal(infile,'spec_width0')[0])))
    sys.stdout.write('   sw(P):  %f ppm\n' % (sw/sfrq))

    sys.stdout.write('Inirect dimension f3:\n')
    dfrq=float(vpar.GetOmegaVal(infile,'f3_freq')[0])
    sw1=float(vpar.GetOmegaVal(infile,'sw1')[0])
    sys.stdout.write('   dfrq:   %f MHz\n' % (dfrq))
    sys.stdout.write('   ni:     %i \n' % (int(nb/2)))
    sys.stdout.write('   at:     %f ms\n' % (int(nb/2)*1/sw1*1000.))
    sys.stdout.write('   sw1:    %f Hz\n' % (sw1))
    sys.stdout.write('   sw1(P): %f ppm\n' % (sw1/dfrq))

    f1180_flg=int(vpar.GetOmegaVal(infile,'f1180_flg')[0])

    #NEED TO: get spectrometer frequency, nuclei and carrier
    nproc=2  #number of dimensions to fourier transform

    temp=float(vpar.GetOmegaVal(infile,'temp')[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(temp))
    vpar.PipeParse('omega',infile.split('.par')[0]+'.bin',ndim,nproc,(2*np,nb),(np,nb/2),('Complex','Complex'),(sw,sw1),(sfrq,dfrq),(vpar.WaterPPM(temp),c13shift,0.0),('1H','13C'),f1180_flg)

    os.system('csh fid.test.com')




def ConvertSpec_nhsqc_omeg(infile,ni=0):

    na=int(vpar.GetOmegaVal(infile,'na')[0])
    nb=int(vpar.GetOmegaVal(infile,'nb')[0])
    nc=int(vpar.GetOmegaVal(infile,'nc')[0])

    np=int(vpar.GetOmegaVal(infile,'block_size')[0])


    ndim=2
    sys.stdout.write('%i\n' % (ndim)) 
    sys.stdout.write('Direct dimension f1:\n')
    sfrq=float(vpar.GetOmegaVal(infile,'dim0_freq')[0])
    sw=float(vpar.GetOmegaVal(infile,'spec_width0')[0])
    sys.stdout.write('   sfrq:   %f MHz\n' % (sfrq))
    sys.stdout.write('   nt:     %i \n' % (na))
    sys.stdout.write('   at:     %f ms\n' % (float(vpar.GetOmegaVal(infile,'dwell_time0')[0])*1000.))
    sys.stdout.write('   np:     %i\n' % (np))
    sys.stdout.write('   sw:     %f Hz\n' % (float(vpar.GetOmegaVal(infile,'spec_width0')[0])))
    sys.stdout.write('   sw(P):  %f ppm\n' % (sw/sfrq))


    if(float(sfrq)>600 and float(sfrq)<610): #on the omega600, nitrogen is on f3 for this seq
        sys.stdout.write('Inirect dimension f3:\n')
        dfrq=float(vpar.GetOmegaVal(infile,'f3_freq')[0])
    else:
        sys.stdout.write('Inirect dimension f2:\n')
        dfrq=float(vpar.GetOmegaVal(infile,'f2_freq')[0])

    sw1=1/(float(vpar.GetOmegaVal(infile,'t1dw')[0])*1E-6)
    sys.stdout.write('   dfrq:   %f MHz\n' % (dfrq))
    sys.stdout.write('   ni:     %i \n' % (int(nc)))
    sys.stdout.write('   at:     %f ms\n' % (int(nc)*1/sw1*1000.))
    sys.stdout.write('   sw1:    %f Hz\n' % (sw1))
    sys.stdout.write('   sw1(P): %f ppm\n' % (sw1/dfrq))

    #f1180_flg=int(vpar.GetOmegaVal(infile,'f1180_flg')[0])
    f1180_flg=1

    #NEED TO: get spectrometer frequency, nuclei and carrier
    nproc=2  #number of dimensions to fourier transform

    print vpar.GetOmegaVal(infile,'temp')

    temp=float(vpar.GetOmegaVal(infile,'temp')[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(temp))
    vpar.PipeParse('omega',infile.split('.par')[0]+'.bin',ndim,nproc,(2*np,2*nc),(np,nc),('Complex','Rance-Kay'),(sw,sw1/2),(sfrq,dfrq),(vpar.WaterPPM(temp),n15shift,0.0),('1H','15N'),f1180_flg)

    os.system('csh fid.test.com')




def ConvertSpec_nhsqc_omeg(infile,ni=0):

    na=int(vpar.GetOmegaVal(infile,'na')[0])
    nb=int(vpar.GetOmegaVal(infile,'nb')[0])
    nc=int(vpar.GetOmegaVal(infile,'nc')[0])

    np=int(vpar.GetOmegaVal(infile,'block_size')[0])


    ndim=2
    sys.stdout.write('%i\n' % (ndim)) 
    sys.stdout.write('Direct dimension f1:\n')
    sfrq=float(vpar.GetOmegaVal(infile,'dim0_freq')[0])
    sw=float(vpar.GetOmegaVal(infile,'spec_width0')[0])
    sys.stdout.write('   sfrq:   %f MHz\n' % (sfrq))
    sys.stdout.write('   nt:     %i \n' % (na))
    sys.stdout.write('   at:     %f ms\n' % (float(vpar.GetOmegaVal(infile,'dwell_time0')[0])*1000.))
    sys.stdout.write('   np:     %i\n' % (np))
    sys.stdout.write('   sw:     %f Hz\n' % (float(vpar.GetOmegaVal(infile,'spec_width0')[0])))
    sys.stdout.write('   sw(P):  %f ppm\n' % (sw/sfrq))


    if(float(sfrq)>600 and float(sfrq)<610): #on the omega600, nitrogen is on f3 for this seq
        sys.stdout.write('Inirect dimension f3:\n')
        dfrq=float(vpar.GetOmegaVal(infile,'f3_freq')[0])
    else:
        sys.stdout.write('Inirect dimension f2:\n')
        dfrq=float(vpar.GetOmegaVal(infile,'f2_freq')[0])

    sw1=1/(float(vpar.GetOmegaVal(infile,'t1dw')[0])*1E-6)
    sys.stdout.write('   dfrq:   %f MHz\n' % (dfrq))
    sys.stdout.write('   ni:     %i \n' % (int(nc)))
    sys.stdout.write('   at:     %f ms\n' % (int(nc)*1/sw1*1000.))
    sys.stdout.write('   sw1:    %f Hz\n' % (sw1))
    sys.stdout.write('   sw1(P): %f ppm\n' % (sw1/dfrq))

    #f1180_flg=int(vpar.GetOmegaVal(infile,'f1180_flg')[0])
    f1180_flg=1

    #NEED TO: get spectrometer frequency, nuclei and carrier
    nproc=2  #number of dimensions to fourier transform

    print vpar.GetOmegaVal(infile,'temp')

    temp=float(vpar.GetOmegaVal(infile,'temp')[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(temp))
    vpar.PipeParse('omega',infile.split('.par')[0]+'.bin',ndim,nproc,(2*np,2*nc),(np,nc),('Complex','Rance-Kay'),(sw,sw1/2),(sfrq,dfrq),(vpar.WaterPPM(temp),n15shift,0.0),('1H','15N'),f1180_flg)

    os.system('csh fid.test.com')


def ConvertSpec_nhsqc_bruk(ni=0):

    na=int(vpar.GetBrukVal('acqus','NS')[0])
    nb=int(vpar.GetBrukVal('acqu2s','TD')[0])
#    nc=int(vpar.GetBrukVal('acqu3s','TD')[0])

    np=int(vpar.GetBrukVal('acqus','TD')[0])

    ndim=2
    sys.stdout.write('%i\n' % (ndim)) 
    sys.stdout.write('Direct dimension f1:\n')
    sfrq=float(vpar.GetBrukVal('acqus','SFO1')[0])
    sw=float(vpar.GetBrukVal('acqus','SW_h')[0])
    sys.stdout.write('   sfrq:   %f MHz\n' % (sfrq))
    sys.stdout.write('   nt:     %i \n' % (na))
    sys.stdout.write('   at:     %f ms\n' % (np/sw))
    sys.stdout.write('   np:     %i\n' % (np))
    sys.stdout.write('   sw:     %f Hz\n' % (sw))
    sys.stdout.write('   sw(P):  %f ppm\n' % (sw/sfrq))

    sys.stdout.write('Inirect dimension f2:\n')
    
    
    dfrq=float(vpar.GetBrukVal('acqu2s','SFO1')[0])
    sw1=float(vpar.GetBrukVal('acqu2s','SW_h')[0])

    sys.stdout.write('   dfrq:   %f MHz\n' % (dfrq))
    sys.stdout.write('   ni:     %i \n' % (int(nb)/2))
    sys.stdout.write('   at:     %f ms\n' % (int(nb)*1/sw1*1000.))
    sys.stdout.write('   sw1:    %f Hz\n' % (sw1))
    sys.stdout.write('   sw1(P): %f ppm\n' % (sw1/dfrq))

    #f1180_flg=int(vpar.GetOmegaVal(infile,'f1180_flg')[0])
    f1180_flg=1

    #NEED TO: get spectrometer frequency, nuclei and carrier
    nproc=2  #number of dimensions to fourier transform

    #print vpar.GetOmegaVal(infile,'temp')

    temp=float(vpar.GetBrukVal('acqu','TEMP')[0])
    print 'water:',vpar.WaterPPM(temp)

    #inflate the total number of points to make it even out at a multiple of 1024
    npAdj=BrukFidAdjust(np)



    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(temp))

    vpar.PipeParse('bruk','ser',ndim,nproc,(npAdj,nb),(np/2,nb/2),('Complex','Rance-Kay'),(sw,sw1),(sfrq,dfrq),(vpar.WaterPPM(temp),n15shift,0.0),('1H','15N'),f1180_flg)

    os.system('csh fid.test.com')


#inflate the total number of points to make it even out at a multiple of 1024
def BrukFidAdjust(np):
    factor=math.ceil(np/2.*4. / 1024.)
    npAdj=1024*factor/2.
    return npAdj


def ConvertSpec_HtoCssdisp(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])


    nz=len(vpar.GetParVarian('./procpar','n',('','ncyc_cp')))
    if(nz==1):
        nz=len(vpar.GetParVarian('./procpar','n',('','time_T2')))

    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(1)
    elif(array[1]=='phase'):
        mode=str(0)
    elif(array[2]=='phase'):
        mode=str(0)

    print 'ni:',ni,'np:',np,'nz:',nz,'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

    print 're-shuffle.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid'
    os.system('re-shuffle.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')
    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final.0 fid.0')
    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final.1 fid.1')
    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final.2 fid.2')
    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final.3 fid.3')

    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Complex','Real'),(sw,sw1,0.0),(sfrq,dfrq,0.0),(vpar.WaterPPM(temp),c13shift,0.0),('1H','13C','ncyc_cp'),f1180_flg,loop=4)    

    os.system('csh fid.test.com')



def ConvertSpec_Forbidden(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    forbid_flg=vpar.GetParVarian('./procpar','n',('','forbid_flg'))
    time_T2=vpar.GetParVarian('./procpar','n',('','time_T2'))

    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(2)
    elif(array[1]=='phase'):
        mode=str(1)
    elif(array[0]=='phase'):
        mode=str(2)

    nz=len(time_T2)
    print 'ni:',ni,'np:',np,'nz:',nz,'flags:',len(forbid_flg),'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

    os.system('RelaxFix2.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+str(len(forbid_flg))+' '+mode+' fid.final fid')

    tof=float(vpar.GetParVarian('./procpar','n',('','tof'))[0])
    tof_me=float(vpar.GetParVarian('./procpar','n',('','tof_me'))[0])

    direct=(tof_me-tof)/sfrq+vpar.WaterPPM(temp)

    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz*len(forbid_flg)),(np/2,ni,nz*len(forbid_flg)),('Complex','Complex','Real'),(sw,sw1,0.0),(sfrq,dfrq,0.0),(direct,c13shift,0.0),('1H','13C','ncyc_cp'),f1180_flg)    

    os.system('csh fid.test.com')



def ConvertSpec_IzIzz(ni=0):

    if(ni==0):
        ni=int(vpar.GetParVarian('./procpar','n',('','ni'))[0])
    np=int(vpar.GetParVarian('./procpar','n',('','np'))[0])
    Iz_flg=vpar.GetParVarian('./procpar','n',('','Iz_flg'))
    time_T1=vpar.GetParVarian('./procpar','n',('','time_T1'))

    array=vpar.GetParVarian('./procpar','n',('','array'))[0].split('"')[1].split(',')
    if(array[0]=='phase'):
        mode=str(1)
    elif(array[2]=='phase'):
        mode=str(0)

    nz=len(time_T1)*len(Iz_flg)
    print 'ni:',ni,'np:',np,'nz:',len(time_T1),'flags:',len(Iz_flg),'mode:',mode

    sw=float(vpar.GetParVarian('./procpar','n',('','sw'))[0])
    sw1=float(vpar.GetParVarian('./procpar','n',('','sw1'))[0])

    sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])
    dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])
    dfrq2=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])
    f1180_flg=vpar.GetParVarian('./procpar','n',('','f1180'))[0]

    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    print 'water:',vpar.WaterPPM(temp)

    c13shift,n15shift=vpar.Cshift(sfrq,dfrq,dfrq2,vpar.WaterPPM(temp))

    os.system('RelaxFix.out '+str(np)+' '+str(ni)+' '+str(nz)+' '+mode+' fid.final fid')

#    tof=float(vpar.GetParVarian('./procpar','n',('','tof'))[0])
#    tof_me=float(vpar.GetParVarian('./procpar','n',('','tof_me'))[0])
#    direct=(tof_me-tof)/sfrq+vpar.WaterPPM(temp)
    direct=vpar.WaterPPM(temp)

    vpar.PipeParse('var','fid.final',3,3,(np,ni*2,nz),(np/2,ni,nz),('Complex','Complex','Real'),(sw,sw1,0.0),(sfrq,dfrq,0.0),(direct,c13shift,0.0),('1H','13C','ncyc_cp'),f1180_flg)    

    os.system('csh fid.test.com')




def ConvertSpec(ni=0,ppmMin='*',ppmMax='*'):
    
    type=vpar.GetSpectrometerType()
    
    if(type=='var'):
        seqfil=vpar.GetParVarian('./procpar','n',('','seqfil'))[0].split('"')[1]
        print seqfil



        if(seqfil=='HtoC_CH3_exchange_600_lek_ILV'):
            ConvertSpec_HtoC(ni=ni)
        elif(seqfil=='CH3_forbiddenDQ_allowed_600_lek'):
            ConvertSpec_Forbidden(ni=ni)
        elif(seqfil=='HtoC_CH3_exchange_600_DC_dfh_v2' or seqfil=='HtoC_CH3_exchange_600_DC_dfh_v2_forAB'):
            ConvertSpec_HtoCssdisp(ni=ni)                         
        elif(seqfil=='N15_CPMG_Rex_NH_fm_600_v6'):
            ConvertSpec_15Ncpmg(ni=ni)
        elif(seqfil=='N15_CPMG_Rex_NH_fm_500'):
            ConvertSpec_15Ncpmg(ni=ni)
        elif(seqfil=='N15_CPMG_Rex_NH_fm_500_v5'):
            ConvertSpec_15Ncpmg(ni=ni)
        elif(seqfil=='N15_CPMG_Rex_NH_fm_500_v6'):
            ConvertSpec_15Ncpmg(ni=ni)

        elif(seqfil=='CT_N_hsqc_LED_lek_600_v2'):
            ConvertSpec_15N_LED(ni=ni)

        elif(seqfil=='hmqc_c13_600_methyl_diffusion_lek'):
            ConvertSpec_13C_diff(ni=ni)


	elif(seqfil=='hsqc_gd_sl_seduce_NHT2_600'):
	    Convertspec_15NHT2(ni=ni)
            
        #Typical order parameter experiments
	elif(seqfil=='N15NOE_lek_pfg_sel_enh_600'):
	    Convertspec_15NHNOE(ni=ni)
	elif(seqfil=='N15T2_lek_pfg_sel_enh_600'):
	    Convertspec_15NHT2b(ni=ni)
	elif(seqfil=='N15T1_lek_pfg_sel_enh_600'):
	    Convertspec_15NHT1(ni=ni)

	elif(seqfil=='N15T1_lek_pfg_sel_enh_600'):
	    Convertspec_15NHT1(ni=ni)

	elif(seqfil=='NHT1_unenhanced_lek_600'):
	    Convertspec_15NHT1unenhanced(ni=ni)

        elif(seqfil=='CH3_1HT2s_600_lek'):
            ConvertSpec_ProtonT2s(ni=ni)
        elif(seqfil=='CH3_T1Z_T1ZZ_lek_600'):
            ConvertSpec_IzIzz(ni=ni)
        elif(seqfil=='hmqc_c13_600_methyl_lek'):
            ConvertSpec_chmqc(ni=ni,ppmMin=ppmMin,ppmMax=ppmMax)
        elif(seqfil=='hsqc_gd_sl_seduce_500'):
            ConvertSpec_nhsqc(ni=ni,ppmMin=ppmMin,ppmMax=ppmMax)
        elif(seqfil=='hsqc_gd_sl_seduce_600'):
            ConvertSpec_nhsqc(ni=ni,ppmMin=ppmMin,ppmMax=ppmMax)
        elif(seqfil=='CT_hsqc_600'):
            ConvertSpec_ct_chsqc(ni=ni,ppmMin=ppmMin,ppmMax=ppmMax)
        elif(seqfil=='N15_CPMG_Rex_NH_trosy_antitrosy_lek_600_v4'):
            ConvertSpec_15NcpmgTrosyAntiTrosy(ni=ni,ppmMin=ppmMin,ppmMax=ppmMax)

        elif(seqfil=='water_sLED_fm_v2_600'):
            ConvertSpec_waterSLED(ni=ni,ppmMin=ppmMin,ppmMax=ppmMax)
        elif(seqfil=='tnnoesy_ajb'):
            ConvertSpec_noesy(ni=ni)

        else:
            print 'Cannot find ',seqfil,' in varian processing library'
            sys.exit(100)
        sys.exit(100)
    if(type=='omega'):
        
        parfile=vpar.GetOmegaParFile()
        test=vpar.GetParOmega(parfile,'n',('','seq_source',))[0].split('/')
        seqfil=test[len(test)-1]

        print 'Sequence file:',seqfil

        if(seqfil=='ajb_hmqc_c13_methyl.s'):
            ConvertSpec_chmqc_omeg(parfile,ni=ni)
        elif(seqfil=='ajb_HtoC_CH3_exchange_600_ILV.s'):
            ConvertSpec_HtoC_omeg(parfile,ni=ni)
        elif(seqfil=='ajb_HtoC_CH3_exchange_750_ILV.s'):
            ConvertSpec_HtoC_omeg(parfile,ni=ni)
        elif(seqfil=='ajb_HtoC_CH3_exchange_950_ILV.s'):
            ConvertSpec_HtoC_omeg(parfile,ni=ni)

        elif(seqfil=='GrEnhHsqc_wfb.s'):
            ConvertSpec_nhsqc_omeg(parfile,ni=ni)
        elif(seqfil=='ajb_GrEnhHsqc_wfb.s'):
            ConvertSpec_nhsqc2_omeg(parfile,ni=ni)

        elif(seqfil=='ajb_N15_CPMG_NH_500.s'):
            ConvertSpec_15Ncpmg_omeg(parfile,ni=ni)
        elif(seqfil=='ajb_N15_CPMG_NH_600.s'):
            ConvertSpec_15Ncpmg_omeg(parfile,ni=ni)
        elif(seqfil=='ajb_N15_CPMG_NH_900.s'):
            ConvertSpec_15Ncpmg_omeg(parfile,ni=ni)
        else:
            print 'Cannot find ',seqfil,' in varian processing library'
            sys.exit(100)
    else:
        seqfil=vpar.GetParBruk('acqu','n',('','PULPROG',))[0].split('<')[1].split('>')[0]
        print seqfil        

        if(seqfil=='15N_CPMG_N_Rex_cw.ajb'):
            ConvertSpec_15Ncpmg_bruk('acqu',ni=ni)
        if(seqfil=='r2_disp_exp_3D.mk'):
            ConvertSpec_15Ncpmg_bruk('acqu',ni=ni)

        elif(seqfil=='hsqcetfpf3gpsi'):
            ConvertSpec_nhsqc_bruk(ni=ni)
        else:
            print 'Cannot find ',seqfil,' in bruker processing library'
            sys.exit(100)
