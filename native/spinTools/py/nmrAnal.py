#!/usr/bin/python

####################################################################
# Functions to quantitatively analyse NMR intensities 
# Function to use depends on the experiment and desired analysis
#
# A.Baldwin 15th Feb 2013
####################################################################

import os,string,math,vpar,sys,numpy,scipy
from baldwinStd import readfile,PathExists


#first two sequences need to be updated on a pulse-sequence by pulse-sequence basis.
#bit of a pain, but this is the only way to get it right I think.
def GetCPMGparams():


    
    type=vpar.GetSpectrometerType()
    if(type=='var'):
        seqfil=vpar.GetParVarian('./procpar','n',('','seqfil'))[0].split('"')[1]
        time_T2=float(vpar.GetParVarian('./procpar','n',('','time_T2'))[0])
        temp=(vpar.GetParVarian('./procpar','n',('','temp'))[0])        


        if(seqfil=='N15_CPMG_Rex_NH_trosy_antitrosy_lek_600_v4'):        
            field=str(int(round(float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0]))))
            pwx_cp=float(vpar.GetParVarian('./procpar','n',('','pwn_cp'))[0])
            taub=float(vpar.GetParVarian('./procpar','n',('','taub'))[0])
            tau_eq=float(vpar.GetParVarian('./procpar','n',('','time_equil'))[0])

            seqfil_cpmg='Trosy_CPMG'
            basis='TrATr_13;estimate_etaz'

#            seqfil_cpmg='PE_CPMG'
#            basis='IphAph_13'

            sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])        
            dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])        

            c13shift,xcar=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(float(temp)))
            nucleus='N15'
            coupled='H1'

            local=[]
            local.append(('Omega',xcar))
            local.append(('JIS',-92))            
            local.append(('DeltaJ',0))            
            local.append(('R0Tr_%s' % field,5.0))            
            local.append(('R0ATr_%s' % field,5.0))            
            local.append(('R1aph_%s' % field,2.0))            
            local.append(('R1iph_%s' % field,1.0))            
            local.append(('DeltaO',1.0))            
            fix=[]
            fix.append('Omega')
            fix.append('JIS')
            fix.append('DeltaJ')
            fix.append('R0ATr_%s' % field)
            fix.append('R1iph_%s' % field)
            fix.append('R1aph_%s' % field)

        if(seqfil=='HtoC_CH3_exchange_600_lek_ILV'):        
            field=str(int(round(float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0]))))
            pwx_cp=float(vpar.GetParVarian('./procpar','n',('','pwc_cp'))[0])
            taub=float(vpar.GetParVarian('./procpar','n',('','taub'))[0])
            #tau_eq=float(vpar.GetParVarian('./procpar','n',('','time_equil'))[0])
            tau_eq=0.000

            seqfil_cpmg='PE_CPMG'
            basis='IphAph_13'

            sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])        
            dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq'))[0])        

            xcar,n15shift=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(float(temp)))
            nucleus='C13'
            coupled='H1'

            local=[]
            local.append(('Omega',xcar))
            local.append(('JIS',125))            
            local.append(('DeltaJ',0))            
            local.append(('R0iph_%s' % field,5.0))            
#            local.append(('R0ATr_%s' % field,5.0))            
            local.append(('R1aph_%s' % field,2.0))            
            local.append(('R1iph_%s' % field,1.0))            
            local.append(('DeltaO',1.0))            
            fix=[]
            fix.append('Omega')
            fix.append('JIS')
            fix.append('DeltaJ')
            fix.append('R0_%s' % field)
            fix.append('R1iph_%s' % field)
            fix.append('R1aph_%s' % field)


        if(seqfil=='N15_CPMG_Rex_NH_fm_500' or seqfil=='N15_CPMG_Rex_NH_fm_500_v6'):        
            field=str(int(round(float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0]))))
            pwx_cp=float(vpar.GetParVarian('./procpar','n',('','pwn_cp'))[0])
            taub=float(vpar.GetParVarian('./procpar','n',('','taub'))[0])
            #tau_eq=float(vpar.GetParVarian('./procpar','n',('','time_equil'))[0])
            tau_eq=0.000

            seqfil_cpmg='PE_CPMG'
            basis='IphAph_13'

            sfrq=float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])        
            dfrq=float(vpar.GetParVarian('./procpar','n',('','dfrq2'))[0])        

            xcar,n15shift=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(float(temp)))
            nucleus='N15'
            coupled='H1'

            local=[]
            local.append(('Omega',xcar))
            local.append(('JIS',125))            
            local.append(('DeltaJ',0))            
            local.append(('R0iph_%s' % field,5.0))            
#            local.append(('R0ATr_%s' % field,5.0))            
            local.append(('R1aph_%s' % field,2.0))            
            local.append(('R1iph_%s' % field,1.0))            
            local.append(('DeltaO',1.0))            
            fix=[]
            fix.append('Omega')
            fix.append('JIS')
            fix.append('DeltaJ')
            fix.append('R0_%s' % field)
            fix.append('R1iph_%s' % field)
            fix.append('R1aph_%s' % field)






    elif(type=='omega'):
        parfile=vpar.GetOmegaParFile()
        test=vpar.GetParOmega(parfile,'n',('','seq_source',))[0].split('/')
        seqfil=test[len(test)-1]
#        parfile=vpar.GetParFile()
        time_T2=float(vpar.GetParOmega(parfile,'n',('','time_T2'))[0])*1E-6

    elif(type=='bruk'):
        seqfil=vpar.GetParBruk('acqu','n',('','PULPROG',))[0].split('<')[1].split('>')[0]
        time_T2=float(vpar.GetParBruk('acqus','n',('','D21'))[0])*2.
        temp=(vpar.GetParBruk('acqu','n',('','TEMP'))[0])
        if(seqfil=='15N_CPMG_N_Rex_cw.ajb'):
            field=str(int(round(float(vpar.GetParBruk('acqus','n',('','BF1'))[0]))))

            pwx_cp=float(vpar.GetParBruk('acqus','n',('','P25'))[0])
            tau_eq=float(vpar.GetParBruk('acqus','n',('','D28'))[0])
            taub=float(vpar.GetParBruk('acqus','n',('','D25'))[0])
            seqfil_cpmg='Trosy_CPMG'
            basis='TrATr_13;estimate_etaz'

#            seqfil_cpmg='PE_CPMG'
#            basis='IphAph_13'

            sfrq=float(vpar.GetBrukVal('acqus','SFO1')[0])
            dfrq=float(vpar.GetBrukVal('acqu2','SFO1')[0])
            c13shift,xcar=vpar.Cshift(sfrq,dfrq,dfrq,vpar.WaterPPM(float(temp)))
            nucleus='N15'
            coupled='H1'

            local=[]
            local.append(('Omega',xcar))
            local.append(('JIS',-92))            
            local.append(('DeltaJ',0))            
            local.append(('R0Tr_%s' % field,5.0))            
            local.append(('R0ATr_%s' % field,5.0))            
            local.append(('R1aph_%s' % field,2.0))            
            local.append(('R1iph_%s' % field,1.0))            
            local.append(('DeltaO',1.0))            
            fix=[]
            fix.append('Omega')
            fix.append('JIS')
            fix.append('DeltaJ')
            fix.append('R0ATr_%s' % field)
            fix.append('R1iph_%s' % field)
            fix.append('R1aph_%s' % field)

            
    #each dataset requires its own local parameters
    loc=[]
    outy=open('catia/ParamSet_'+field+'_'+temp+'.inp','w')
    loc.append('catia/ParamSet_'+field+'_'+temp+'.inp')
    outy.write('format = (')
    for i in range(len(local)):
        if(i!=0):
            outy.write(';')
        outy.write('%s' % local[i][0])
    outy.write(')\n')
    outy.write('* = (')
    for i in range(len(local)):
        if(i!=0):
            outy.write(';')
        outy.write('%s' % str(local[i][1]))        
    outy.write(')\n')
    outy.close()
    
    outy=open('catia/ParamGlobal.inp','w')
    outy.write('kex=1000.\n')
    outy.write('pb=0.02\n')
    outy.close()
    glob=[]
    glob.append('catia/ParamGlobal.inp')

    return loc,glob,fix,field,temp,pwx_cp,taub,time_T2,tau_eq,xcar,seqfil_cpmg,basis,nucleus,coupled


def WriteCatiaFile(datasets,loc,glob,fix):
    outy=open('catia/OneFieldFit.catia','w')
    
    for i in range(len(datasets)):
        outy.write('ReadDataset(%s)  #Data summary file \n' % datasets[i])
    for i in range(len(loc)):
        outy.write('ReadParam_Local(%s)      #Local parameters initial\n' % loc[i])
    for i in range(len(glob)):
        outy.write('ReadParam_Global(%s)   #Global parameters initial\n' % glob[i])

    outy.write('ReadParam(Omega;test.ft2.list;0;1)   #peak list\n')
    if(os.path.exists('./catia/DeltaOmega.inp')):
        outy.write('ReadParam(DeltaO;./catia/DeltaOmega.inp;0;1)#delta omegas\n')

    outy.write('# Fix all the static parameters\n')

    for i in range(len(fix)):
        outy.write('FreeLocalParam(all;%s;false)\n' % fix[i])

    outy.write('FreeLocalParam(all;DeltaO;false)\n')


#    outy.write('SetGlobalParam(kex;250)\n')
#    outy.write('SetGlobalParam(pb;0.05)\n')

    PathExists(('./catia/OutPut',))
    outy.write('# Deal with global parameters\n')
    outy.write('FreeGlobalParam(kex;true)\n')
    outy.write('FreeGlobalParam(pb;true)\n')
    outy.write('# Minimize\n')
    outy.write('echo(\\n)\n')
    outy.write('#Minimize()\n')
    outy.write('Minimize(print=y;tol=1E-3;maxiter=100)\n')
    outy.write('#\n')
    outy.write('#  // Print some files for plotting\n')
    outy.write('PrintParam(./catia/OutPut/GlobalParam.fit;global)\n')
    outy.write('PrintParam(./catia/OutPut/DeltaOmega.fit;DeltaO)\n')
    outy.write('PrintData(./catia/OutPut/)\n')
    outy.write('echo(\n)\n')
    outy.write('ChiSq(all;all)\n')
    outy.write('exit(0)\n')
    outy.close()


def WriteCatiaDataset(am,var):
    PathExists(('catia','catia/dataset'))    

    loc,glob,fix,field,temp,pwx_cp,taub,time_T2,tau_eq,xcar,seqfil_cpmg,basis,nuc,coupled=GetCPMGparams()

    datasets=[]



    if(var=='n'):
        outy=open('./catia/dataset/Dataset_'+field+'_'+temp+'_Trosy.inp','w')
        datasets.append('./catia/dataset/Dataset_'+field+'_'+temp+'_Trosy.inp')
    else:
        outy=open('./catia/dataset/Dataset_'+field+'_'+temp+'_Trosy_temp.inp','w')
        datasets.append('./catia/dataset/Dataset_'+field+'_'+temp+'_Trosy_temp.inp')

    outy.write('ID=15N CPMG Trosy @ %s and %s\n' % (field,temp))
    outy.write('sfrq = %s\n' % (field))

    if(var=='y'):
        #outy.write('severalTemperatures = true \n')
        #outy.write('multipleTemperatures = yes \n')
        #outy.write('kexType = arrhenius \n')
        #outy.write('pbType = arrhenius \n')
        #outy.write('rateType = arrhenius \n')
        outy.write('deltaOmegaType = linear \n')
#        outy.write('deltaOmegaType = standard \n')


    outy.write('temperature = %s\n' % (temp))
    outy.write('nucleus = %s\n' % nuc)
    outy.write('couplednucleus = %s\n' % coupled)
    outy.write('pwx_cp = %sE-6\n' % pwx_cp)#in seconds
    outy.write('taub = %s \n' % taub)        #in seconds
    outy.write('time_T2 = %s \n' % time_T2)  #in seconds
    outy.write('time_equil = %s\n' % tau_eq)#in seconds
    outy.write('xcar = %.2f\n' % xcar)

    outy.write('minerror = (1.%;0.3/s)\n')

    outy.write('seqfil = %s\n' % seqfil_cpmg)
    outy.write('basis = (%s)\n' % basis)

#    outy.write('seqfil = PE_CPMG\n')
#    outy.write('basis = (IphAph_13)\n')

    outy.write('format = (0;1;2)\n')
    outy.write('DataDirectory = %s \n' % ('./fuda/')   )
    outy.write('Data = (\n')
    for i in range(len(am)):
        outy.write('[%s;%s];\n' % (am[i],am[i]+'.out.cpmg'))
    outy.write(')\n')
    outy.close()
    
    return datasets,loc,glob,fix


def CPMGAnal():

    type=vpar.GetSpectrometerType()
    if(type=='var'):
        TIME_T2=float(vpar.GetParVarian('./procpar','n',('','time_T2'))[0])
    elif(type=='omega'):
        parfile=vpar.GetParFile()
        TIME_T2=float(vpar.GetParOmega(parfile,'n',('','time_T2'))[0])*1E-6
    elif(type=='bruk'):
        seqfil=vpar.GetParBruk('acqu','n',('','PULPROG',))[0].split('<')[1].split('>')[0]
        if(seqfil=='15N_CPMG_N_Rex_cw.ajb'):
            TIME_T2=float(vpar.GetParBruk('acqus','n',('','D21'))[0])*2.
        if(seqfil=='r2_disp_exp_3D.mk'):
            TIME_T2=30E-3



    else:
        print 'No rule for spectrometer type',type

    vals=[]
    filey=os.listdir('fuda')
    for i in range(len(filey)):
        if(len(filey[i].split('.out'))>1):
            if(len(filey[i].split('.cpmg'))==1):
                vals.append('fuda/'+filey[i])
            
    gnu=open("cpmg.gnu",'w')
    gnu.write("set encoding iso_8859_1\n")
    gnu.write("set xlabel \'{/Symbol n}_{CPMG} (Hz)\' \n")
    gnu.write("set ylabel \'R@_2^{eff} (s^{-1})\' \n")
    gnu.write("set xrange[*:*]\n")
    gnu.write("set yrange[*:*]\n")
    gnu.write("set term post enh mono solid 18\n")
    gnu.write("set out \'%s\' \n" % ("cpmg.ps"))

    gnun=open("cpmgDisp.gnu",'w')
    gnun.write("set encoding iso_8859_1\n")
    gnun.write("set xlabel \'{/Symbol n}_{CPMG} (Hz)\' \n")
    gnun.write("set ylabel \'R@_2^{eff} (s^{-1})\' \n")
    gnun.write("set xrange[*:*]\n")
    gnun.write("set yrange[*:*]\n")
    gnun.write("set term post enh mono solid 18\n")
    gnun.write("set out \'%s\' \n" % ("cpmgDisp.ps"))

    xmin=0.
    xmax=1.
    ymin=100.
    ymax=1.

    yminDisp=100.
    ymaxDisp=1.

    am=[]
    for val in vals: #for each peak datafile
        # Get name of the peak
        Name=string.split(string.strip(val),'.')[0]
        Name=Name.replace("fuda/","")
        am.append(Name)
        # Load data file
        inputfile = open(val,'r')
        lines=inputfile.readlines()
        Data=[]
        dublicate=[]
        for line in lines:
            if ( len(string.split(line)) > 2 ):
                if string.split(line)[1]=='f01(ppm)':
                    Noffset=float(string.split(line)[2])
                if string.split(line)[1]=='f02(ppm)':
                    Hoffset=float(string.split(line)[2])
            #
            if not ( line[0] == "#"):
                #
                # Read the data.
                temp=string.split(line)
                if ( math.fabs(float(temp[0]))<1e-6):
                    Ref=[0,float(temp[1]),float(temp[2])]
                else:
                    Data.append([float(temp[0]),float(temp[1]),float(temp[2])])
                    
                for d in range(len(Data)-1):
                    if ( math.fabs(float(temp[0])-Data[d][0]) < 1e-6 ):
                        dublicate.append(math.pow(float(temp[1])-Data[d][1],2.))
        #
        # Get dublicate data (Flemming's original spelling mistake)
        if(len(dublicate)==0):
            StdErr=0.3
        else:
            StdErr=0.
            for d in range(len(dublicate)):
                StdErr+=dublicate[d]
            StdErr=math.sqrt(StdErr)/len(dublicate)
        #
        # Calculate the R2 and field
        ofn=open(val+".cpmg",'w')
        ofn.write("#%11s%15s%13s\n" % ('nu_cpmg(Hz)','R2(1/s)','Esd(R2)'))
        
        for j in range(len(Data)):
            d=Data[j]
            field=d[0]/TIME_T2

            R=math.log(math.fabs(Ref[1]/d[1]))/TIME_T2
            Esd=math.sqrt( math.pow(StdErr/Ref[1],2.)+math.pow(StdErr/d[1],2.))/TIME_T2


            if(j==0):
                resXmin=field
                resXmax=field
                resYmax=R
                resYmin=R
                resYmaxerr=Esd
                resYminerr=Esd


            if(field<xmin):
                xmin=field
            if(field>xmax):
                xmax=field
            if(R+Esd>ymax):
                if(R+Esd<100):
                    ymax=R+Esd
            if(R-Esd<ymin):
                if(R-Esd>-100):
                    ymin=R-Esd
                    
            if(field<resXmin):
                resXmin=field
                resYmin=R
                resYmaxerr=Esd
            if(field>resXmax):
                resXmax=field
                resYmax=R
                resYmaxerr=Esd
            


            ofn.write(" %11.4e%15.6e%13.6e\n" % (field,R,Esd))
        gnu.write("plot \'%s\' u 1:2:3 t \'%s \@ (%.1fppm; %.1fppm)\' w e pt 6 ps 1.5 \n" % (val+".cpmg",Name,Noffset,Hoffset))

        #print resXmax,resYmin,resYmax,resYmin
        if( (resYmax+resYmaxerr)< (resYmin-resYminerr) ):
            print 'Dispersion!',val
            if(yminDisp>resYmin-resYminerr):
                yminDisp=resYmin-resYminerr
            if(ymaxDisp<resYmax-resYmaxerr):
                ymaxDisp=resYmax+resYmaxerr
            gnun.write("plot \'%s\' u 1:2:3 t \'%s \@ (%.1fppm; %.1fppm)\' w e pt 6 ps 1.5 \n" % (val+".cpmg",Name,Noffset,Hoffset))
        ofn.close()
        inputfile.close()

    gnu.close()
    gnun.close()

    gnu=open('cpmg.gnu')#adjust x and y ranges
    gnu2=open('cpmg2.gnu','w')
    for line in gnu.readlines():
        if(len(line.split('xrange'))>1):
            gnu2.write('set xrange[%f:%f]\n' % (xmin,xmax*1.05))
        elif(len(line.split('yrange'))>1):
            gnu2.write('set yrange[%f:%f]\n' % (ymin,ymax))
        else:
            gnu2.write(line)
    gnu2.close()
    gnu.close()

    gnu=open('cpmgDisp.gnu')#adjust x and y ranges
    gnu2=open('cpmgDisp2.gnu','w')
    for line in gnu.readlines():
        if(len(line.split('xrange'))>1):
            gnu2.write('set xrange[%f:%f]\n' % (xmin,xmax*1.05))
        elif(len(line.split('yrange'))>1):
            gnu2.write('set yrange[%f:%f]\n' % (yminDisp,ymaxDisp))
        else:
            gnu2.write(line)
    gnu2.close()
    gnu.close()

    os.system('gnuplot cpmg2.gnu')
    os.system('sed -e "/gnulinewidth 5.000 def/s//gnulinewidth 12.000 def/" cpmg.ps > _t')
    os.system('mv _t cpmg.ps')

    os.system('gnuplot cpmgDisp2.gnu')
    os.system('sed -e "/gnulinewidth 5.000 def/s//gnulinewidth 12.000 def/" cpmgDisp.ps > _t')
    os.system('mv _t cpmgDisp.ps')


    #setup Catia files

    datasets,loc,glob,fix=WriteCatiaDataset(am,'n')
    WriteCatiaFile(datasets,loc,glob,fix)
    
    if(os.uname()[0]=='Darwin'):
        os.system('catia_Darwin_i386 < catia/OneFieldFit.catia')
    else:
        os.system('catia_Linux_x86_64 < catia/OneFieldFit.catia')


    MakeGnuCPMG('catia/OutPut')





def GetParamLocal(infile,param):
    array=readfile(infile)
    for i in range(len(array)):
        if(len(array[i])>2):
            if(array[i][0]=='#'):
                if(array[i][1]==param):
                    return array[i][2],array[i][3]
    return 0,0

def GetDataParam(infile,param):
    array=readfile(infile)
    outy=[]
    for i in range(len(array)):
        if(len(array[i])>1):
            if(array[i][0]=='#'+param):
                outy.append(array[i][1])
    return outy

def GetParamGlobal(infile,param):
    array=readfile(infile)
    for i in range(len(array)):
        if(len(array[i])>2):
            if(array[i][0]==param):
                return array[i][1],array[i][2]
    return 0,0


def MakeLabel(outdir,gnun,val,temps):
    if(len(temps)==0):
        deltaO,deltaOerr=GetParamLocal(outdir+'/'+val,'DeltaO')
        kex,kexerr=GetParamGlobal(outdir+'/GlobalParam.fit','kex')
        pb,pberr=GetParamGlobal(outdir+'/GlobalParam.fit','pb')
        

        if(kexerr=='fixed'):
            gnun.write('set label sprintf(\'k_{ex} (s^{-1}): %s +/- fixed \',%s) at graph 0.05,0.15\n' % ('%.2f',kex))
        else:
            gnun.write('set label sprintf(\'k_{ex} (s^{-1}): %s +/- %s \',%s,%s) at graph 0.05,0.15\n' % ('%.2f','%.2f',kex,kexerr))

        if(pberr=='fixed'):
            gnun.write('set label sprintf(\'P_b : %s +/- fixed \',%s) at graph 0.05,0.1\n' % ('%.2f',float(pb)*100))
        else:
            gnun.write('set label sprintf(\'P_b : %s +/- %s \',%s,%s) at graph 0.05,0.1\n' % ('%.2f','%.2f',float(pb)*100,float(pberr)*100))

        if(deltaOerr=='fixed'):
            gnun.write('set label sprintf(\'{/Symbol Dw}(ppm): %s +/- fixed \',%s) at graph 0.05,0.05\n' % ('%.2f',deltaO))
        else:
            gnun.write('set label sprintf(\'{/Symbol Dw}(ppm): %s +/- %s \',%s,%s) at graph 0.05,0.05\n' % ('%.2f','%.2f',deltaO,deltaOerr))

    else:
        deltaOa,deltaOaerr=GetParamLocal(outdir+'/'+val,'DeltaO_a')
        deltaOb,deltaOberr=GetParamLocal(outdir+'/'+val,'DeltaO_b')
        dHb,dHberr=GetParamGlobal(outdir+'/GlobalParam.fit','deltaHb')
        dSb,dSberr=GetParamGlobal(outdir+'/GlobalParam.fit','deltaSb')
        dHab,dHaberr=GetParamGlobal(outdir+'/GlobalParam.fit','deltaHab')
        dSab,dSaberr=GetParamGlobal(outdir+'/GlobalParam.fit','deltaSab')

        gnun.write('set label sprintf(\'{/Symbol Dw}(ppm)\') at graph 0.25,0.95\n')
        gnun.write('set label sprintf(\'P_b\') at graph 0.5,0.95\n')
        gnun.write('set label sprintf(\'k_{ex} (s^{-1})\') at graph 0.7,0.95\n')
        for i in range(len(temps)):
            
            tempK=float(temps[i])+273.19

            PLANCKS_H=    6.62606896E-34   # kg m^2 / s 
            BOLTZMANN_KB= 1.3806504E-23    # kg m^2 / K s^2 
            MOLAR_GAS_R=  8.314472         # kg m^2 / K mol s^2           
            A = 3000 * tempK;
            RT = MOLAR_GAS_R * tempK;
            
            
            dG=(float(dHb)-float(dSb)*tempK)
                #Keq=math.exp(-dG/RT)
            dGab=float(dHab)-tempK*float(dSab)
            dGba=dGab-dG

            kab=A*math.exp(-dGab/RT)
            kba=A*math.exp(-dGba/RT)
                
            kex=kab+kba
            pb=kba/kex
            
            gnun.write('set label sprintf(\'%s^oC:\') at graph 0.05,%f\n' % (temps[i],0.9-i*0.05))
            gnun.write('set label sprintf(\'%s\',%f) at graph 0.25,%f\n' % ('%.2f',float(deltaOb)+float(deltaOa)*(float(temps[i])),0.9-i*0.05))
            gnun.write('set label sprintf(\'%s\',%f) at graph 0.5,%f\n' % ('%.2f',pb*100,0.9-i*0.05))
            gnun.write('set label sprintf(\'%s\',%f) at graph 0.7,%f\n' % ('%.2f',kex,0.9-i*0.05))


#print out a summary of the experimental dispersions
def MakeGnuCPMG(outdir,temps=[]):
    tag=outdir.split('/')[1]
    print 'Making pretty outputs for',tag
    gnun=open("norm.gnu",'w')
    gnun.write("set encoding iso_8859_1\n")
    gnun.write("set xlabel \'{/Symbol n}_{CPMG} (Hz)\' \n")
    gnun.write("set ylabel \'R@_2^{eff} (s^{-1})\' \n")
    gnun.write("set term post enh color solid 18\n")
    gnun.write("set size square\n")
    gnun.write("set key outside\n")
    gnun.write("set yrange[-1:14]\n")
    gnun.write("set xrange[*:*]\n")
    gnun.write("set out \'%s\' \n" % ("fitnorm."+tag+".ps"))

    gnu=open("plot.gnu",'w')
    gnu.write("set encoding iso_8859_1\n")
    gnu.write("set xlabel \'{/Symbol n}_{CPMG} (Hz)\' \n")
    gnu.write("set ylabel \'R@_2^{eff} (s^{-1})\' \n")
    gnu.write("set size square\n")
    gnu.write("set key outside\n")
    gnu.write("set yrange[-1:14]\n")
    gnu.write("set xrange[*:*]\n")
    gnu.write("set term post enh color solid 18\n")
    gnu.write("set out \'%s\' \n" % ("fitcpmg."+tag+".ps"))

    files=os.listdir("./"+outdir+"/")
    vals=[]
    for i in range(len(files)):
        if(len(files[i].split('.dat'))>1):
            vals.append(files[i])

    #initialise limits for plotting
    minny1=0.
    maxxy1=1.
    minny0=0.
    maxxy0=1.
    xmax=1.
    xmin=0.
    for k in range(len(vals)):
        # Get name of the peak
        Name=outdir+"/"+string.split(string.strip(vals[k]),'.')[0]

        #Produce files whose last point is subtracted
        inputfile = open(outdir+"/"+vals[k],'r')
        lines=inputfile.readlines()

        tick=0 #analyse the data to break into blocks
        Data=[]
        dataset=[]
        for line in lines:
            if not ( line[0] == "#"):
                temp=string.split(line)
                if(len(temp)!=0):
                    dataset.append(temp)
                else:
                    tick=1
            if(tick==1):
                if(len(dataset)>0):
                    Data.append(dataset)
                    dataset=[]
                tick=0

        outfile=open(Name+".norm",'w') #subtract last number    
        for j in range(len(Data)):
            for i in range(len(Data[j])):
                x1=float(Data[j][i][0])
                y1=float(Data[j][i][1])-float(Data[j][len(Data[j])-1][3])
                s1=float(Data[j][i][2])
                y2=float(Data[j][i][3])-float(Data[j][len(Data[j])-1][3])
                
                if(y1+s1 > maxxy1):
                    maxxy1=y1+s1
                if(y1-s1 < minny1):
                    minny1=y1-s1
                if(float(Data[j][i][1])+s1 > maxxy0):
                    maxxy0=float(Data[j][i][1])+s1
                if(float(Data[j][i][1])-s1 < minny0):
                    minny0=float(Data[j][i][1])-s1
                if(x1<xmin):
                    xmin=x1
                if(x1>xmax):
                    xmax=x1

                outfile.write('%f\t%f\t%f\t%f\n' % (x1,y1,s1,y2))
            outfile.write('\n\n')
        outfile.close()        


        field=GetDataParam(outdir+'/'+vals[k],'Field:')
        temp=GetDataParam(outdir+'/'+vals[k],'Temperature:')

        MakeLabel(outdir,gnun,vals[k],temps)

        #produce plotting script
        Name=Name.replace(outdir+'/',"")
        gnun.write("plot ")
        for i in range(len(Data)):
            if(i!=0):
                gnun.write(',')
            gnun.write("\'%s\' i %i u 1:2:3 t \'%s^oC \' w e pt 6 ps 2 lt %i, '' i %i u 1:4 t \'\' w lines lt %i" % (outdir+"/"+Name+".norm",i,Name+"/"+field[i]+" "+temp[i],i+1,i,i+1))
        gnun.write('\n')
        gnun.write('unset label\n')

        MakeLabel(outdir,gnu,vals[k],temps)

        gnu.write("plot ")
        for i in range(len(Data)):
            if(i!=0):
                gnu.write(',')
            gnu.write("\'%s\' i %i u 1:2:3 t \'%s^oC \' w e pt 6 ps 2 lt %i, '' i %i u 1:4 t \'\' w lines lt %i" % (outdir+"/"+Name+".dat",i,Name+"/"+field[i]+" "+temp[i],i+1,i,i+1))
        gnu.write('\n')
        gnu.write('unset label\n')

    gnun.close()
    gnu.close()

    test=open('norm.gnu','r')
    outy=open('norm2.gnu','w')
    for line in test.readlines():
        if(len(line.split('yrange'))>1):
            outy.write('set yrange[%f:%f]\n' % (minny1,maxxy1))
        elif(len(line.split('xrange'))>1):
            outy.write('set xrange[%f:%f]\n' % (xmin,xmax*1.05))
        else:
            outy.write(line)
    outy.close()

    test=open('plot.gnu','r')
    outy=open('plot2.gnu','w')
    for line in test.readlines():
        if(len(line.split('yrange'))>1):
            outy.write('set yrange[%f:%f]\n' % (minny0,maxxy0))
        elif(len(line.split('xrange'))>1):
            outy.write('set xrange[%f:%f]\n' % (xmin,xmax*1.05))
        else:
            outy.write(line)
    outy.close()


    os.system("gnuplot norm2.gnu")
    os.system("gnuplot plot2.gnu")

    os.system('sed -e "/gnulinewidth 5.000 def/s//gnulinewidth 12.000 def/" fitnorm.'+tag+'.ps > _t')
    os.system('/bin/mv _t fitnorm.'+tag+'.ps')

    os.system('sed -e "/gnulinewidth 5.000 def/s//gnulinewidth 12.000 def/" fitcpmg.'+tag+'.ps > _t')
    os.system('/bin/mv _t fitcpmg.'+tag+'.ps')
    return temp[0]






# Take input file, and take ratios between lines i and i+1
# for analysis of forbidden experiment
# 
# Tugarinov et al, JACS 129,1743-1750 2007
#
# by A. Baldwin
# November 2009
###################################################################
def ForbidAnal():

    #first setup gnuplot file that will do the fitting (old skool!)
    hbar    =  1.05457148E-34 #m^2 kg s-1
    RHH     =  1.81E-10       #m
    P2cosR6 = -0.5/RHH**3     # m^-3
    gammaH  =  2.67512896E8
    mu0d4pi =  1E-7
    CONN    =  hbar**2 * gammaH**4 * P2cosR6**2 *9/10 * mu0d4pi **2/1E9  #constant, to make units ns


    array=vpar.GetParVarian('./procpar','n',('','array'))
    forbid_flg=vpar.GetParVarian('./procpar','n',('','forbid_flg'))
    time_T2=vpar.GetParVarian('./procpar','n',('','time_T2'))

    tag=0
    for test in (array[0].split(',')):
        if(test=='"time_T2"'):
            tag=1
        if(test=='"forbid_flg"' and tag!=1):
            print 'Experiment is not time_T2,forbid_flg'
            sys.exit(100)
    
    if(len(forbid_flg)!=2):
        print 'Problem! Need two forbid_flg elements'
        sys.exit(100)

    if(forbid_flg[0]!='"y"'):
        print 'Problem! Need forbid_flg to go y,n'
        sys.exit(100)
            

    xmax=max(time_T2)

    #from the paper:
    #Ia/Ib= -0.5 nu * tanh (sqrt((nu^2+delta^2)*T))/(sqrt(nu^2+delta^2)-delta*tanh(sqrt((nu^2+delta^2)*T)))
    # where:
    # nu    = 9/10 (P2cos^theta/R3)^2 Saxis^2tauC * gammaH^4 * hbar^2
    # delta = -4/9 Rext 

    gnu=open('gnu.gp','w')
    gnu.write('set term post eps enh color solid 20\n')
    gnu.write('set output \'forbid.eps\'\n')
    gnu.write('\n')
    gnu.write('set fit errorvariables\n')
    gnu.write('hbar   = 1.05E-34\n')
    gnu.write('RHH    = 1.813E-10\n')
    gnu.write('P2cosR6=-0.5/RHH**3\n')
    gnu.write('gammaH =2.67512896E8\n')
    gnu.write('mu0d4pi=1E-7\n')
    gnu.write('CONN= hbar**2 * gammaH**4 * P2cosR6**2 *9/10 * mu0d4pi **2/1E9\n')
    gnu.write('f1(x,nu,de)=0.5 * nu *tanh( ((nu**2+de**2)**0.5)*x)/ ((nu**2+de**2)**0.5 - de*tanh(((nu**2+de**2)**0.5)*x))\n')
    gnu.write('set size square\n')
    gnu.write('set key left\n')
    gnu.write('set xlabel \'T2 (s^{-1})\'\n')
    gnu.write('set ylabel \'|Ia/Ib|\'\n')
    gnu.write('unset key\n')
    gnu.write('set xrange[0:%s]\n' % (float(xmax)*1.05))

    files=os.listdir('fuda')
    filearr=[]
    for file in files:
        test=file.split('.out')
        if(len(test)>1):
            if(len(test[1])==0):
                filearr.append(file)


    #loop over all read in files in fuda directory
    for file in range(len(filearr)):

        print file,filearr[file]   
        name=string.split(filearr[file],'.')[0]
   
        print name
        inputfile = open('fuda/'+filearr[file],'r')
        lines=inputfile.readlines()
        input=[]
        for line in lines:
            if not ( line[0] == "#"):
                # Read the data.
                temp=string.split(line)
                input.append(temp)

        if(len(input)!=len(time_T2)*2):
            print 'Danger: number of time_T2 points does not jive!'
            sys.exit(100)

        data=[]
        for i in range(len(input)/2): #assume list goes forbid/not forbid and time_T2 is correct
            xval=float(time_T2[i])
            yval=(float(input[2*i][1]))/float(input[2*i+1][1])*-1.0/1.0  #adjust the forbid='y' #assumes scans are equally weighted
            yerr=yval*math.sqrt( (float(input[2*i][2])/float(input[2*i][1]))**2  + (float(input[2*i+1][2])/float(input[2*i+1][1]))**2)
            data.append((xval,yval,yerr))
        outputfile = open('fuda/'+filearr[file]+".conv",'w')
        for i in range(len(data)):
            outputfile.write('%e\t%e\t%e\n' % (data[i][0],data[i][1],data[i][2]))
        outputfile.close()

        gnu.write('set title \'%s\'\n' % name) #make individual output eps files
        gnu.write('set output \'fuda/%s.eps\'\n' % name)
        gnu.write('nu1=100;de1=50;\n')
        gnu.write('fit f1(x,nu1,de1) \'fuda/%s.out.conv\' u 1:(abs($2)):3 via nu1,de1\n' % name)
        gnu.write('set label sprintf("{/Symbol h} %.1f +/- %.1f s^{-1} {/Symbol d} %0.1f +/- %0.1f s^{-1}",nu1,nu1_err,de1,de1_err) at graph 0,graph 0.2\n')
        gnu.write('set label sprintf("S^2{/Symbol t}_c %.2f +/- %.2f ns R_{ext} %0.1f +/- %0.1f s^{-1}",nu1/CONN,nu1_err/CONN,de1*(-9/4),de1_err*(9/4)) at graph 0,graph 0.1\n')
        gnu.write('plot \'fuda/%s.out.conv\' u 1:2:3 w err lt 1,f1(x,nu1,de1) lt 1\n' % name)
        gnu.write('unset label\n')
   
    gnu.close()
    os.system('rm fit.log')
    os.system('gnuplot gnu.gp')
    os.system('arraygraph.py 4 6 20 0 20 0 `ls fuda/*.eps`') #make the obligatury summary.pdf


    fitty=readfile('fit.log')
    cnt=0
    tag=0
    list=[]
    for i in range(len(fitty)):
        if(len(fitty[i])>0):
            if(fitty[i][0]=='Final'):
                cnt+=1
                tag=1
        if(tag==1):
            if(len(fitty[i])>0):
                if(fitty[i][0]=='nu1'):
                    nu1  = float(fitty[i][2])/CONN*1.0
                    nu1e = float(fitty[i][4])/CONN*1.0
                if(fitty[i][0]=='de1'):
                    de1  = float(fitty[i][2])*-9.0/4.0
                    de1e = float(fitty[i][4])*9.0/4.0
                    list.append((cnt,nu1,nu1e,de1,de1e))
                    tag=0
         
#if (len(list)!=len(sys.argv)-1):
#   print 'Problem - lists do not match up'
#   sys.exit(100)


    forbb=open('forbid.out','w') #make a nice file with all the values in
    forbb.write('#\tS^2tau_c (ns) \t\t\tRext (s-1)\n')
    ass=[]

    for file in range(len(filearr)):
        name=string.split(filearr[file],'.')[0]
        ass.append(((string.split(name,'C-H')[0]),name,list[file][1],list[file][2],list[file][3],list[file][4]))

    ass=sorted(ass,key=lambda ass: ass[0])
    for i in range(len(ass)):
        forbb.write('%s\t%e\t%e\t%e\t%e\n' % (ass[i][1],ass[i][2],ass[i][3],ass[i][4],ass[i][5]))
    forbb.close()




###################################################################

def IzIzzAnal():

    array=vpar.GetParVarian('./procpar','n',('','array'))
    Iz_flg=vpar.GetParVarian('./procpar','n',('','Iz_flg'))
    time_T1=vpar.GetParVarian('./procpar','n',('','time_T1'))

    tag=0
    for test in (array[0].split(',')):
        if(test=='"time_T1"'):
            tag=1
        if(test=='"Iz_flg"' and tag!=1):
            print 'Experiment is not time_T1,Iz_flg. Need to code the exception!'
            sys.exit(100)
    
    if(len(Iz_flg)!=2):
        print 'Problem! Need two forbid_flg elements. Code exception!'
        sys.exit(100)

    if(Iz_flg[0]!='"y"'):
        print 'Problem! Need forbid_flg to go y,n. Code exception!'
        sys.exit(100)
            

    xmax=max(time_T1)

    files=os.listdir('fuda')
    filearr=[]
    for file in files:
        test=file.split('.out')
        if(len(test)>1):
            if(len(test[1])==0):
                filearr.append(file)


    sfrq=int( math.ceil(float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])))
    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])



    gnu=open('gnu.gp','w')
    gnu.write('set term post eps enh color solid \'Arial\'20\n')
    gnu.write('set fit errorvariables\n')
    gnu.write('f(x,m,r)=m*exp(-r*x)\n')
    gnu.write('set xrange[0:%f]\n' % (float(xmax)*1.05*1000.))
    gnu.write('set size square\n')
    gnu.write('set xlabel \'T (ms)\'\n')
    gnu.write('set ylabel \'I/I_0\'\n')

    for i in range(len(filearr)):
        
        array=readfile('fuda/'+filearr[i])
        outy=open('fuda/'+filearr[i]+'.conv','w')
        data=[]
        for j in range(len(array)):
            if(array[j][0]!='#'):
                if(len(array[j])==3):
                    data.append((array[j][0],array[j][1],array[j][2]))
        if(len(data)!=len(time_T1)*2):
            print 'Data shapes do not match'
            print len(time_T1),len(data)

            sys.exit(100)
        for j in range(len(data)/2):
            ex1=time_T1[j]
            ey1d=data[2*j][1]
            ey1e=data[2*j][2]
            ey2d=data[2*j+1][1]
            ey2e=data[2*j+1][2]
            outy.write('%s\t%s\t%s\t%s\t%s\n' % (ex1,ey1d,ey1e,ey2d,ey2e))
        outy.close()


        name=filearr[i].split('.')[0]
        
        array=readfile('fuda/'+filearr[i])
        intEst=1E6
        for j in range(len(array)):
            if(len(array[j])==3):
                try:
                    intEst=float(array[j][1])
                    break
                except:
                    pass

        gnu.write('set title \'CarbonIzIzz %s Field: %s Temp: %s\'\n' % (name,str(sfrq),str(temp))) #make individual output eps files
        gnu.write('set output \'fuda/%s.eps\'\n' % name)

        gnu.write('m=%f;r=1E-1\n' % (intEst))
        gnu.write('fit f(x,m,r) \'fuda/%s.out.conv\' u 1:(abs($2)) via m,r\n' % (name))
        gnu.write('set label sprintf("%s Iz %s +/- %s s^{-1}",r,r_err) font "Arial,12" at graph 0.05,graph 0.2\n' % (name,'%.1f','%.1f'))

        gnu.write('n=%f;s=1E-1\n' % (intEst))
        gnu.write('fit f(x,n,s) \'fuda/%s.out.conv\' u 1:(abs($4)) via n,s\n' % (name))
        gnu.write('set label sprintf("%s Izz %s +/- %s s^{-1}",s,s_err) font "Arial,12" at graph 0.05,graph 0.15\n' % (name,'%.1f','%.1f'))

        gnu.write('plot \'fuda/%s.out.conv\' u ($1*1000):($2/m):($3/m) noti w err lt 1,(f(x/1000,m,r)/m) ti \'Iz\' lt 1,\\\n' % name)
        gnu.write('\'fuda/%s.out.conv\' u ($1*1000):($4/n):($5/n) noti w err lt 2,(f(x/1000,n,s)/n) ti \'Izz\' lt 2\n' % name)
        gnu.write('unset label\n')

        

    gnu.close()


    os.system('rm fit.log')
    os.system('gnuplot gnu.gp')

    #os.system('arraygraph.py 4 6 20 0 20 0 `ls fuda/*.eps`') #make the obligatury summary.pdf


    fitty=readfile('fit.log')
    cnt=0
    tag=0
    list=[]
    for i in range(len(fitty)):
        if(len(fitty[i])>0):
            if(fitty[i][0]=='Final'):
                cnt+=1
                tag=1
        if(tag==1):
            if(len(fitty[i])>4):
                if(fitty[i][0]=='m'):
                    m  = float(fitty[i][2])
                    me = float(fitty[i][4])
                if(fitty[i][0]=='r'):
                    r  = float(fitty[i][2])
                    re = float(fitty[i][4])
                if(fitty[i][0]=='n'):
                    n  = float(fitty[i][2])
                    ne = float(fitty[i][4])
                if(fitty[i][0]=='s'):
                    s  = float(fitty[i][2])
                    se = float(fitty[i][4])
                    list.append((cnt,m,me,r,re,n,ne,s,se))
                    tag=0

         
    forbb=open('IzIzz.out','w') #make a nice file with all the values in
    forbb.write('#\tI01 \t\t\tRIz (s-1)\t\tI02  \t\t RIzz (s-1)\n')
    ass=[]
    for file in range(len(filearr)):
        name=string.split(filearr[file],'.')[0]
        ass.append(((string.split(name,'C-H')[0]),name,list[file][1],list[file][2],list[file][3],list[file][4],list[file][5],list[file][6],list[file][7],list[file][8]))

    ass=sorted(ass,key=lambda ass: ass[0])
    for i in range(len(ass)):
        forbb.write('%s\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n' % (ass[i][1],ass[i][2],ass[i][3],ass[i][4],ass[i][5],ass[i][6],ass[i][7],ass[i][8],ass[i][9]))
    forbb.close()






###################################################################
# ProtonT2s analysis
def ProtonT2sAnal():

    files=os.listdir('fuda')
    filearr=[]
    for file in files:
        test=file.split('.out')
        if(len(test)>1):
            if(len(test[1])==0):
                filearr.append(file)


    sfrq=int( math.ceil(float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])))
    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    time_T2=vpar.GetParVarian('./procpar','n',('','time_T2'))
    xmax=max(time_T2)

    gnu=open('gnu.gp','w')
    gnu.write('set term post eps enh color solid \'Arial\'20\n')
    gnu.write('set fit errorvariables\n')
    gnu.write('f(x,m,r)=m*exp(-r*x)\n')
    gnu.write('set xrange[0:%f]\n' % (float(xmax)*1.05*1000.))
    gnu.write('set size square\n')
    gnu.write('unset key\n')
    gnu.write('set xlabel \'T (ms)\'\n')
    gnu.write('set ylabel \'I/I_0\'\n')

    for i in range(len(filearr)):
        name=filearr[i].split('.')[0]
        
        array=readfile('fuda/'+filearr[i])
        intEst=1E6
        for j in range(len(array)):
            if(len(array)==3):
                try:
                    intEst=float(array[i][1])
                except:
                    pass

        gnu.write('set title \'ProtonR2 %s Field: %s Temp: %s\'\n' % (name,str(sfrq),str(temp))) #make individual output eps files
        gnu.write('set output \'fuda/%s.eps\'\n' % name)

        gnu.write('m=%f;r=1E2\n' % (intEst))
        gnu.write('fit f(x,m,r) \'fuda/%s.out\' u 1:(abs($2)):3 via m,r\n' % (name))
        gnu.write('set label sprintf("%s proton T_2 %s +/- %s s^{-1}",r,r_err) font "Arial,12" at graph 0.05,graph 0.2\n' % (name,'%.1f','%.1f'))
        gnu.write('plot \'fuda/%s.out\' u ($1*1000):($2/m):($3/m) w err lt 1,(f(x/1000,m,r)/m) lt 1\n' % name)
        gnu.write('unset label\n')
    gnu.close()

    os.system('rm fit.log')
    os.system('gnuplot gnu.gp')
    os.system('arraygraph.py 4 6 20 0 20 0 `ls fuda/*.eps`') #make the obligatury summary.pdf


    fitty=readfile('fit.log')
    cnt=0
    tag=0
    list=[]
    for i in range(len(fitty)):
        if(len(fitty[i])>0):
            if(fitty[i][0]=='Final'):
                cnt+=1
                tag=1
        if(tag==1):
            if(len(fitty[i])>0):
                if(fitty[i][0]=='m'):
                    nu1  = float(fitty[i][2])
                    nu1e = float(fitty[i][4])
                if(fitty[i][0]=='r'):
                    de1  = float(fitty[i][2])
                    de1e = float(fitty[i][4])
                    list.append((cnt,nu1,nu1e,de1,de1e))
                    tag=0
         
    forbb=open('protonT2s.out','w') #make a nice file with all the values in
    forbb.write('#\tI0 \t\t\tR (s-1)\n')
    ass=[]
    for file in range(len(filearr)):
        name=string.split(filearr[file],'.')[0]
        ass.append(((string.split(name,'C-H')[0]),name,list[file][1],list[file][2],list[file][3],list[file][4]))

    ass=sorted(ass,key=lambda ass: ass[0])
    for i in range(len(ass)):
        forbb.write('%s\t%e\t%e\t%e\t%e\n' % (ass[i][1],ass[i][2],ass[i][3],ass[i][4],ass[i][5]))
    forbb.close()


def ConvFloat(array):
    arr=[]
    for i in range(len(array)):
        arr.append(float(array[i]))
    return arr

###################################################################
# ProtonT2s analysis
def OrdNT2():

    files=os.listdir('fuda')
    filearr=[]
    for file in files:
        test=file.split('.out')
        if(len(test)>1):
            if(len(test[1])==0):
                filearr.append(file)


    sfrq=int( math.ceil(float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])))
    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])

    ncyc=vpar.GetParVarian('./procpar','n',('','ncyc'))
    ncyc=numpy.array(ConvFloat(ncyc))
    pwn=float(vpar.GetParVarian('./procpar','n',('','pwn'))[0])
    time_T2=ncyc*(32.0*pwn*1E-6 + 32.0*450.0e-6)
    xmax=max(time_T2)



    gnu=open('gnu.gp','w')
    gnu.write('set term post eps enh color solid 20\n')
    gnu.write('set fit errorvariables\n')
    gnu.write('f(x,m,r)=m*exp(-r*x)\n')
    gnu.write('set xrange[0:%f]\n' % (float(xmax)*1.05*1000.))
    gnu.write('set size square\n')
    gnu.write('unset key\n')
    gnu.write('set xlabel \'T (ms)\'\n')
    gnu.write('set ylabel \'I/I_0\'\n')

    for i in range(len(filearr)):
        name=filearr[i].split('.')[0]
        
        array=readfile('fuda/'+filearr[i])
        intEst=1E6
        for j in range(len(array)):
            if(len(array)==3):
                try:
                    intEst=float(array[i][1])
                except:
                    pass

        gnu.write('set title \'NitrogenR2 %s Field: %s Temp: %s\'\n' % (name,str(sfrq),str(temp))) #make individual output eps files
        gnu.write('set output \'fuda/%s.eps\'\n' % name)

        gnu.write('m=%f;r=1E2\n' % (intEst))
        gnu.write('fit f(x,m,r) \'fuda/%s.out\' u 1:2:3 via m,r\n' % (name))
        gnu.write('set label sprintf("%s nitrogen R_2 %s +/- %s s^{-1}",r,r_err) font "Arial,12" at graph 0.05,graph 0.2\n' % (name,'%.1f','%.1f'))
        gnu.write('plot \'fuda/%s.out\' u ($1*1000):($2/m):($3/m) w err lt 1,(f(x/1000,m,r)/m) lt 1\n' % name)
        gnu.write('unset label\n')
    gnu.close()

    os.system('rm fit.log')
    os.system('gnuplot gnu.gp')
    os.system('arraygraph.py 4 6 20 0 20 0 `ls fuda/*.eps`') #make the obligatury summary.pdf


    fitty=readfile('fit.log')
    cnt=0
    tag=0
    list=[]
    for i in range(len(fitty)):
        if(len(fitty[i])>0):
            if(fitty[i][0]=='Final'):
                cnt+=1
                tag=1
        if(tag==1):
            if(len(fitty[i])>0):
                if(fitty[i][0]=='m'):
                    nu1  = float(fitty[i][2])
                    nu1e = float(fitty[i][4])
                if(fitty[i][0]=='r'):
                    de1  = float(fitty[i][2])
                    de1e = float(fitty[i][4])
                    list.append((cnt,nu1,nu1e,de1,de1e))
                    tag=0
         
    forbb=open('nitrogenR2.out','w') #make a nice file with all the values in
    forbb.write('#\tI0 \t\t\tR (s-1)\n')
    ass=[]
    for file in range(len(filearr)):
        name=string.split(filearr[file],'.')[0]
        ass.append(((string.split(name,'C-H')[0]),name,list[file][1],list[file][2],list[file][3],list[file][4]))

    ass=sorted(ass,key=lambda ass: ass[0])
    for i in range(len(ass)):
        forbb.write('%s\t%e\t%e\t%e\t%e\n' % (ass[i][1],ass[i][2],ass[i][3],ass[i][4],ass[i][5]))
    forbb.close()


###################################################################
# ProtonT1s analysis
def OrdNT1():

    files=os.listdir('fuda')
    filearr=[]
    for file in files:
        test=file.split('.out')
        if(len(test)>1):
            if(len(test[1])==0):
                filearr.append(file)


    sfrq=int( math.ceil(float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])))
    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    ncyc=vpar.GetParVarian('./procpar','n',('','ncyc'))
    ncyc=numpy.array(ConvFloat(ncyc))
    pw_shpss=float(vpar.GetParVarian('./procpar','n',('','pw_shpss'))[0])
    time_T2 = ncyc*(pw_shpss*1E-6 + 2.0*2.5e-3);


    xmax=max(time_T2)

    gnu=open('gnu.gp','w')
    gnu.write('set term post eps enh color solid 20\n')
    gnu.write('set fit errorvariables\n')
    gnu.write('f(x,m,r)=m*exp(-r*x)\n')
    gnu.write('set xrange[0:%f]\n' % (float(xmax)*1.05*1000.))
    gnu.write('set size square\n')
    gnu.write('unset key\n')
    gnu.write('set xlabel \'T (ms)\'\n')
    gnu.write('set ylabel \'I/I_0\'\n')

    for i in range(len(filearr)):
        name=filearr[i].split('.')[0]
        
        array=readfile('fuda/'+filearr[i])
        intEst=1E6
        for j in range(len(array)):
            if(len(array)==3):
                try:
                    intEst=float(array[i][1])
                except:
                    pass

        gnu.write('set title \'NitrogenR1 %s Field: %s Temp: %s\'\n' % (name,str(sfrq),str(temp))) #make individual output eps files
        gnu.write('set output \'fuda/%s.eps\'\n' % name)

        gnu.write('m=%f;r=1E2\n' % (intEst))
        gnu.write('fit f(x,m,r) \'fuda/%s.out\' u 1:2:3 via m,r\n' % (name))
        gnu.write('set label sprintf("%s nitrogen R_1 %s +/- %s s^{-1}",r,r_err) font "Arial,12" at graph 0.05,graph 0.2\n' % (name,'%.1f','%.1f'))
        gnu.write('plot \'fuda/%s.out\' u ($1*1000):($2/m):($3/m) w err lt 1,(f(x/1000,m,r)/m) lt 1\n' % name)
        gnu.write('unset label\n')
    gnu.close()

    os.system('rm fit.log')
    os.system('gnuplot gnu.gp')
    os.system('arraygraph.py 4 6 20 0 20 0 `ls fuda/*.eps`') #make the obligatury summary.pdf


    fitty=readfile('fit.log')
    cnt=0
    tag=0
    list=[]
    for i in range(len(fitty)):
        if(len(fitty[i])>0):
            if(fitty[i][0]=='Final'):
                cnt+=1
                tag=1
        if(tag==1):
            if(len(fitty[i])>0):
                if(fitty[i][0]=='m'):
                    nu1  = float(fitty[i][2])
                    nu1e = float(fitty[i][4])
                if(fitty[i][0]=='r'):
                    de1  = float(fitty[i][2])
                    de1e = float(fitty[i][4])
                    list.append((cnt,nu1,nu1e,de1,de1e))
                    tag=0
         
    forbb=open('nitrogenR1.out','w') #make a nice file with all the values in
    forbb.write('#\tI0 \t\t\tR (s-1)\n')
    ass=[]
    for file in range(len(filearr)):
        name=string.split(filearr[file],'.')[0]
        ass.append(((string.split(name,'C-H')[0]),name,list[file][1],list[file][2],list[file][3],list[file][4]))

    ass=sorted(ass,key=lambda ass: ass[0])
    for i in range(len(ass)):
        forbb.write('%s\t%e\t%e\t%e\t%e\n' % (ass[i][1],ass[i][2],ass[i][3],ass[i][4],ass[i][5]))
    forbb.close()



###################################################################
# ProtonT1s analysis
def OrdNT1unenhanced():

    files=os.listdir('fuda')
    filearr=[]
    for file in files:
        test=file.split('.out')
        if(len(test)>1):
            if(len(test[1])==0):
                filearr.append(file)


    sfrq=int( math.ceil(float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])))
    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    ncyc=vpar.GetParVarian('./procpar','n',('','ncyc'))
    ncyc=numpy.array(ConvFloat(ncyc))
    time_T2 = ncyc*(2.0*12.5e-3);


    xmax=max(time_T2)

    gnu=open('gnu.gp','w')
    gnu.write('set term post eps enh color solid 20\n')
    gnu.write('set fit errorvariables\n')
    gnu.write('f(x,m,r)=m*(1-exp(-r*x))\n')
    gnu.write('set xrange[0:%f]\n' % (float(xmax)*1.05*1000.))
    gnu.write('set size square\n')
    gnu.write('unset key\n')
    gnu.write('set xlabel \'T (ms)\'\n')
    gnu.write('set ylabel \'I/I_0\'\n')

    for i in range(len(filearr)):
        name=filearr[i].split('.')[0]
        
        array=readfile('fuda/'+filearr[i])
        intEst=1E6
        for j in range(len(array)):
            if(len(array)==3):
                try:
                    intEst=float(array[i][1])
                except:
                    pass

        gnu.write('set title \'NitrogenR1 %s Field: %s Temp: %s\'\n' % (name,str(sfrq),str(temp))) #make individual output eps files
        gnu.write('set output \'fuda/%s.eps\'\n' % name)

        gnu.write('m=%f;r=1\n' % (intEst))
        gnu.write('fit f(x,m,r) \'fuda/%s.out\' u 1:2:3 via m,r\n' % (name))
        gnu.write('set label sprintf("%s nitrogen R_1 %s +/- %s s^{-1}",r,r_err) font "Arial,12" at graph 0.05,graph 0.2\n' % (name,'%.1f','%.1f'))
        gnu.write('plot \'fuda/%s.out\' u ($1*1000):($2/m):($3/m) w err lt 1,(f(x/1000,m,r)/m) lt 1\n' % name)
        gnu.write('unset label\n')
    gnu.close()

    os.system('rm fit.log')
    os.system('gnuplot gnu.gp')
    os.system('arraygraph.py 4 6 20 0 20 0 `ls fuda/*.eps`') #make the obligatury summary.pdf


    fitty=readfile('fit.log')
    cnt=0
    tag=0
    list=[]
    for i in range(len(fitty)):
        if(len(fitty[i])>0):
            if(fitty[i][0]=='Final'):
                cnt+=1
                tag=1
        if(tag==1):
            if(len(fitty[i])>0):
                if(fitty[i][0]=='m'):
                    nu1  = float(fitty[i][2])
                    nu1e = float(fitty[i][4])
                if(fitty[i][0]=='r'):
                    de1  = float(fitty[i][2])
                    de1e = float(fitty[i][4])
                    list.append((cnt,nu1,nu1e,de1,de1e))
                    tag=0
         
    forbb=open('nitrogenR1.out','w') #make a nice file with all the values in
    forbb.write('#\tI0 \t\t\tR (s-1)\n')
    ass=[]
    for file in range(len(filearr)):
        name=string.split(filearr[file],'.')[0]
        ass.append(((string.split(name,'C-H')[0]),name,list[file][1],list[file][2],list[file][3],list[file][4]))

    ass=sorted(ass,key=lambda ass: ass[0])
    for i in range(len(ass)):
        forbb.write('%s\t%e\t%e\t%e\t%e\n' % (ass[i][1],ass[i][2],ass[i][3],ass[i][4],ass[i][5]))
    forbb.close()



###################################################################
# ProtonT2s analysis
def OrdNNOE():

    files=os.listdir('fuda')
    filearr=[]
    for file in files:
        test=file.split('.out')
        if(len(test)>1):
            if(len(test[1])==0):
                filearr.append(file)


    sfrq=int( math.ceil(float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])))
    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    time_T2=vpar.GetParVarian('./procpar','n',('','ncyc'))
    xmax=max(time_T2)

    gnu=open('gnu.gp','w')
    gnu.write('set term post eps enh color solid \'Arial\'20\n')
    gnu.write('set fit errorvariables\n')
    gnu.write('f(x,m,r)=m*exp(-r*x)\n')
    gnu.write('set xrange[0:%f]\n' % (float(xmax)*1.05*1000.))
    gnu.write('set size square\n')
    gnu.write('unset key\n')
    gnu.write('set xlabel \'T (ms)\'\n')
    gnu.write('set ylabel \'I/I_0\'\n')

    for i in range(len(filearr)):
        name=filearr[i].split('.')[0]
        
        array=readfile('fuda/'+filearr[i])
        intEst=1E6
        for j in range(len(array)):
            if(len(array)==3):
                try:
                    intEst=float(array[i][1])
                except:
                    pass

        gnu.write('set title \'ProtonR2 %s Field: %s Temp: %s\'\n' % (name,str(sfrq),str(temp))) #make individual output eps files
        gnu.write('set output \'fuda/%s.eps\'\n' % name)

        gnu.write('m=%f;r=1E2\n' % (intEst))
        gnu.write('fit f(x,m,r) \'fuda/%s.out\' u 1:(abs($2)):3 via m,r\n' % (name))
        gnu.write('set label sprintf("%s proton T_2 %s +/- %s s^{-1}",r,r_err) font "Arial,12" at graph 0.05,graph 0.2\n' % (name,'%.1f','%.1f'))
        gnu.write('plot \'fuda/%s.out\' u ($1*1000):($2/m):($3/m) w err lt 1,(f(x/1000,m,r)/m) lt 1\n' % name)
        gnu.write('unset label\n')
    gnu.close()

    os.system('rm fit.log')
    os.system('gnuplot gnu.gp')
    os.system('arraygraph.py 4 6 20 0 20 0 `ls fuda/*.eps`') #make the obligatury summary.pdf


    fitty=readfile('fit.log')
    cnt=0
    tag=0
    list=[]
    for i in range(len(fitty)):
        if(len(fitty[i])>0):
            if(fitty[i][0]=='Final'):
                cnt+=1
                tag=1
        if(tag==1):
            if(len(fitty[i])>0):
                if(fitty[i][0]=='m'):
                    nu1  = float(fitty[i][2])
                    nu1e = float(fitty[i][4])
                if(fitty[i][0]=='r'):
                    de1  = float(fitty[i][2])
                    de1e = float(fitty[i][4])
                    list.append((cnt,nu1,nu1e,de1,de1e))
                    tag=0
         
    forbb=open('protonT2s.out','w') #make a nice file with all the values in
    forbb.write('#\tI0 \t\t\tR (s-1)\n')
    ass=[]
    for file in range(len(filearr)):
        name=string.split(filearr[file],'.')[0]
        ass.append(((string.split(name,'C-H')[0]),name,list[file][1],list[file][2],list[file][3],list[file][4]))

    ass=sorted(ass,key=lambda ass: ass[0])
    for i in range(len(ass)):
        forbb.write('%s\t%e\t%e\t%e\t%e\n' % (ass[i][1],ass[i][2],ass[i][3],ass[i][4],ass[i][5]))
    forbb.close()



def FitFuncDiff(x,b):
    Deff=x[0]
    A0=x[1]
    grads=b[0]
    data=b[1]
    b=b[2]

    Ycalc=A0*numpy.exp(-grads**2.*b*Deff)
    print data,A0,A0*numpy.exp(-grads)
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


def DiffAnal():
    ni=int( math.ceil(float(vpar.GetParVarian('./procpar','n',('','ni'))[0])))
    if(ni==1):
        DiffAnal1D()
    else:
        DiffAnal2D()


def DiffAnal1D():

    import diffAnal

    files=os.listdir('fuda')
    filearr=[]
    for file in files:
        test=file.split('.out')
        if(len(test)>1):
            if(len(test[1])==0):
                filearr.append(file)

    import numpy

    sfrq=int( math.ceil(float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])))
    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    gradt=float(vpar.GetParVarian('./procpar','n',('','gt2'))[0])*2.  #gradient duration in seconds
    grads=numpy.array(vpar.GetParVarianFloat('./procpar','n',('','gzlvl2')))
    delta=float(vpar.GetParVarian('./procpar','n',('','Big_delta'))[0])  #gradient duration in seconds
    
    xmax=max(grads)

    #normalise constants
    DAC=0.002     # G cm-1 DAC-1
    yg=2.67513E4  # *rad s-1 G-1
    bfac=(yg*(gradt))**2*(delta-(gradt)/3)

    gnu=open('gnu.gp','w')
    gnu.write('set term post eps enh color solid \'Arial\'20\n')
    gnu.write('set fit errorvariables\n')
    gnu.write('f(x,m,r)=m*exp(-%f*r*x**2)\n' % bfac)
    gnu.write('set xrange[0:%f]\n' % (float(xmax)*1.05*DAC))
    gnu.write('set size square\n')
    gnu.write('unset key\n')
    gnu.write('set xlabel \'G (G cm-1)\'\n')
    gnu.write('set ylabel \'I/I_0\'\n')
    from scipy.optimize import leastsq

    for i in range(len(filearr)):
        name=filearr[i].split('.')[0]
        
        array=readfile('fuda/'+filearr[i])
        intEst=1E6
        data=[]
        for j in range(len(array)):
            print array[j]
            if(len(array[j])==3):
                try:
                    intEst=float(array[j][1])
                    data.append(numpy.abs(float(array[j][1])))
                except:
                    pass

        x0=leastsq(FitFuncDiff,[1E-7,intEst],args=[grads*DAC,numpy.array(data),bfac,])
        #print x0
        #print grads*DAC
        #sys.exit(100)

        
        gnu.write('set title \'Diffusion %s Field: %s Temp: %s\'\n' % (name,str(sfrq),str(temp))) #make individual output eps files
        gnu.write('set output \'fuda/%s.eps\'\n' % name)

        gnu.write('m=%e;r=%e\n' % (x0[0][1],x0[0][0]))
        gnu.write('fit f(x,m,r) \'fuda/%s.out\' u ($1*%f):(abs($2)):3 via m,r\n' % (name,DAC))
        gnu.write('set label sprintf("%s Diffusion %s +/- %s s^{-1}",r,r_err) font "Arial,12" at graph 0.05,graph 0.2\n' % (name,'%.1e','%.1e'))
        gnu.write('set label sprintf("%s Intensity %s +/- %s s^{-1}",m,m_err) font "Arial,12" at graph 0.05,graph 0.15\n' % (name,'%.1f','%.1f'))
        gnu.write('plot \'fuda/%s.out\' u ($1*%f):(abs($2/m)):($3/m) w err lt 1,(f(x,m,r)/m) lt 1\n' % (name,DAC))
        gnu.write('unset label\n')
    gnu.close()

    os.system('rm fit.log')
    os.system('gnuplot gnu.gp')
    os.system('arraygraph.py 4 6 20 0 20 0 `ls fuda/*.eps`') #make the obligatury summary.pdf


    fitty=readfile('fit.log')
    cnt=0
    tag=0
    list=[]
    for i in range(len(fitty)):
        if(len(fitty[i])>0):
            if(fitty[i][0]=='Final'):
                cnt+=1
                tag=1
        if(tag==1):
            if(len(fitty[i])>0):
                if(fitty[i][0]=='m'):
                    nu1  = float(fitty[i][2])
                    nu1e = float(fitty[i][4])
                if(fitty[i][0]=='r'):
                    de1  = float(fitty[i][2])
                    de1e = float(fitty[i][4])
                    list.append((cnt,nu1,nu1e,de1,de1e))
                    tag=0
         
    forbb=open('diffusion.out','w') #make a nice file with all the values in
    forbb.write('#\tI0 \t\t\tD (cm^2s-1)\n')
    ass=[]
    for file in range(len(filearr)):
        name=string.split(filearr[file],'.')[0]
        ass.append(((string.split(name,'C-H')[0]),name,list[file][1],list[file][2],list[file][3],list[file][4]))

    ass=sorted(ass,key=lambda ass: ass[0])
    for i in range(len(ass)):
        forbb.write('%s\t%e\t%e\t%e\t%e\n' % (ass[i][1],ass[i][2],ass[i][3],ass[i][4],ass[i][5]))
    forbb.close()



def DiffAnal2D():

    files=os.listdir('fuda')
    filearr=[]
    for file in files:
        test=file.split('.out')
        if(len(test)>1):
            if(len(test[1])==0):
                filearr.append(file)

    import numpy

    sfrq=int( math.ceil(float(vpar.GetParVarian('./procpar','n',('','sfrq'))[0])))
    temp=float(vpar.GetParVarian('./procpar','n',('','temp'))[0])
    gradt=float(vpar.GetParVarian('./procpar','n',('','gt2'))[0])*2.  #gradient duration in seconds
    grads=numpy.array(vpar.GetParVarianFloat('./procpar','n',('','gzlvl2')))
    delta=float(vpar.GetParVarian('./procpar','n',('','Big_delta'))[0])  #gradient duration in seconds
    
    xmax=max(grads)

    #normalise constants
    DAC=0.002     # G cm-1 DAC-1
    yg=2.67513E4  # *rad s-1 G-1
    bfac=(yg*(gradt))**2*(delta-(gradt)/3)

    gnu=open('gnu.gp','w')
    gnu.write('set term post eps enh color solid \'Arial\'20\n')
    gnu.write('set fit errorvariables\n')
    gnu.write('f(x,m,r)=m*exp(-%f*r*x**2)\n' % bfac)
    gnu.write('set xrange[0:%f]\n' % (float(xmax)*1.05*DAC))
    gnu.write('set size square\n')
    gnu.write('unset key\n')
    gnu.write('set xlabel \'G (G cm-1)\'\n')
    gnu.write('set ylabel \'I/I_0\'\n')
    from scipy.optimize import leastsq

    for i in range(len(filearr)):
        name=filearr[i].split('.')[0]
        
        array=readfile('fuda/'+filearr[i])
        intEst=1E6
        data=[]
        for j in range(len(array)):
            print array[j]
            if(len(array[j])==3):
                try:
                    intEst=float(array[j][1])
                    data.append(numpy.abs(float(array[j][1])))
                except:
                    pass

        x0=leastsq(FitFuncDiff,[1E-7,intEst],args=[grads*DAC,numpy.array(data),bfac,])
        #print x0
        #print grads*DAC
        #sys.exit(100)

        
        gnu.write('set title \'Diffusion %s Field: %s Temp: %s\'\n' % (name,str(sfrq),str(temp))) #make individual output eps files
        gnu.write('set output \'fuda/%s.eps\'\n' % name)

        gnu.write('m=%e;r=%e\n' % (x0[0][1],x0[0][0]))
        gnu.write('fit f(x,m,r) \'fuda/%s.out\' u ($1*%f):(abs($2)):3 via m,r\n' % (name,DAC))
        gnu.write('set label sprintf("%s Diffusion %s +/- %s s^{-1}",r,r_err) font "Arial,12" at graph 0.05,graph 0.2\n' % (name,'%.1e','%.1e'))
        gnu.write('set label sprintf("%s Intensity %s +/- %s s^{-1}",m,m_err) font "Arial,12" at graph 0.05,graph 0.15\n' % (name,'%.1f','%.1f'))
        gnu.write('plot \'fuda/%s.out\' u ($1*%f):(abs($2/m)):($3/m) w err lt 1,(f(x,m,r)/m) lt 1\n' % (name,DAC))
        gnu.write('unset label\n')
    gnu.close()

    os.system('rm fit.log')
    os.system('gnuplot gnu.gp')
    os.system('arraygraph.py 4 6 20 0 20 0 `ls fuda/*.eps`') #make the obligatury summary.pdf


    fitty=readfile('fit.log')
    cnt=0
    tag=0
    list=[]
    for i in range(len(fitty)):
        if(len(fitty[i])>0):
            if(fitty[i][0]=='Final'):
                cnt+=1
                tag=1
        if(tag==1):
            if(len(fitty[i])>0):
                if(fitty[i][0]=='m'):
                    try:
                        nu1  = float(fitty[i][2])
                        nu1e = float(fitty[i][4])
                    except:
                        pass
                if(fitty[i][0]=='r'):
                    try:
                        de1  = float(fitty[i][2])
                        de1e = float(fitty[i][4])
                        list.append((filearr[cnt-1],nu1,nu1e,de1,de1e))
                        print filearr[cnt],nu1,nu1e,de1,de1e
                    except:
                        pass
                    tag=0
     

    forbb=open('diffusion.out','w') #make a nice file with all the values in
    forbb.write('#\tI0 \t\t\tD (cm^2s-1)\n')
    ass=[]
    for file in range(len(list)):
        name=string.split(list[file][0],'.')[0]
        try:
            ass.append(((string.split(name,'C-H')[0]),name,list[file][1],list[file][2],list[file][3],list[file][4]))
        except:
            try:
                ass.append(((string.split(name,'N-H')[0]),name,list[file][1],list[file][2],list[file][3],list[file][4]))
            except:
                 pass
    ass=sorted(ass,key=lambda ass: ass[0])
    for i in range(len(ass)):
        forbb.write('%s\t%e\t%e\t%e\t%e\n' % (ass[i][1],ass[i][2],ass[i][3],ass[i][4],ass[i][5]))
    forbb.close()









def nmrAnal():
    print 'Analysing NMR data:'
    type=vpar.GetSpectrometerType()
    if(type=='var'):
        seqfil=vpar.GetParVarian('./procpar','n',('','seqfil'))[0].split('"')[1]
        #print seqfil
        if(seqfil=='HtoC_CH3_exchange_600_lek_ILV'):
            CPMGAnal()
        elif(seqfil=='CH3_forbiddenDQ_allowed_600_lek'):
            ForbidAnal()
            ConvertSpec_Forbidden(ni=ni)
        elif(seqfil=='HtoC_CH3_exchange_600_DC_dfh_v2'):
            CPMGAnal()
        elif(seqfil=='N15_CPMG_Rex_NH_fm_600_v6'):
            CPMGAnal()
        elif(seqfil=='N15_CPMG_Rex_NH_fm_500'):
            CPMGAnal()
        elif(seqfil=='N15_CPMG_Rex_NH_fm_500_v6'):
            CPMGAnal()
        elif(seqfil=='CT_N_hsqc_LED_lek_600_v2'):
            DiffAnal()
        elif(seqfil=='CH3_1HT2s_600_lek'):
            ProtonT2sAnal()
        elif(seqfil=='CH3_T1Z_T1ZZ_lek_600'):
            IzIzzAnal()
        elif(seqfil=='hmqc_c13_600_methyl_lek'):
            pass
        elif(seqfil=='hsqc_gd_sl_seduce_600'):
            pass
        elif(seqfil=='CT_hsqc_600'):
            pass
        elif(seqfil=='N15_CPMG_Rex_NH_trosy_antitrosy_lek_600_v4'):
            CPMGAnal()

        elif(seqfil=='N15T2_lek_pfg_sel_enh_600'):
            OrdNT2()
        elif(seqfil=='N15T1_lek_pfg_sel_enh_600'):
            OrdNT1()
        elif(seqfil=='NHT1_unenhanced_lek_600'):
            OrdNT1unenhanced()

        elif(seqfil=='N15NOE_lek_pfg_sel_enh_600'):
            OrdNNOE()



        else:
            print 'Cannot find ',seqfil,' in varian processing library'
            sys.exit(100)
    elif(type=='omega'):
        
        parfile=vpar.GetOmegaParFile()
        test=vpar.GetParOmega(parfile,'n',('','seq_source',))[0].split('/')
        seqfil=test[len(test)-1]

        print 'Sequence file:',seqfil

        if(seqfil=='ajb_hmqc_c13_methyl.s'):
            pass
        elif(seqfil=='ajb_HtoC_CH3_exchange_600_ILV.s'):
            CPMGAnal()
        elif(seqfil=='ajb_HtoC_CH3_exchange_750_ILV.s'):
            CPMGAnal()
        elif(seqfil=='ajb_HtoC_CH3_exchange_950_ILV.s'):
            CPMGAnal()
        else:
            print 'Cannot find ',seqfil,' in omega processing library'
            sys.exit(100)

    elif(type=='bruk'):
        seqfil=vpar.GetParBruk('acqu','n',('','PULPROG',))[0].split('<')[1].split('>')[0]
        print seqfil        
        if(seqfil=='15N_CPMG_N_Rex_cw.ajb'):
            CPMGAnal()
        if(seqfil=='r2_disp_exp_3D.mk'):
            CPMGAnal()
        else:
            print 'Cannot find ',seqfil,' in bruker processing library'
            sys.exit(100)
