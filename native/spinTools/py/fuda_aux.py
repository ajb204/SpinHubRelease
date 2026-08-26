#!/usr/bin/python

#####################################################################
# Auxilary functions to setup fuda input scripts
#
# A.Baldwin 15th Feb 2013

import vpar,sys

def MakeFudaInput(peaklist,radiusH,radiusC,fac,outfile='test.ft2',fudafile='param.fuda',zcoor='ncyc_cp'):
    print 'bullshit'
    delayFac=1.0
    type=vpar.GetSpectrometerType()

    if(type=='var'):
        seqfil=vpar.GetParVarian('./procpar','n',('','seqfil'))[0].split('"')[1]
        print seqfil

        if(seqfil=='HtoC_CH3_exchange_600_lek_ILV'):
            zcoor='ncyc_cp'
#        elif(seqfil=='CH3_forbiddenDQ_allowed_600_lek'):
#            ConvertSpec_Forbidden(ni=ni)
        elif(seqfil=='HtoC_CH3_exchange_600_DC_dfh_v2' or seqfil=='HtoC_CH3_exchange_600_DC_dfh_v2_forAB'):
            zcoor='ncyc_cp'
        elif(seqfil=='hmqc_c13_600_methyl_diffusion_lek'):
            zcoor='gzlvl5'
        elif(seqfil=='N15_CPMG_Rex_NH_fm_600_v6'):
            zcoor='ncyc'
        elif(seqfil=='N15_CPMG_Rex_NH_fm_500'):
            zcoor='ncyc'
        elif(seqfil=='N15_CPMG_Rex_NH_fm_500_v6'):
            zcoor='ncyc'
        elif(seqfil=='CT_N_hsqc_LED_lek_600_v2'):
            zcoor='gzlvl2'
        elif(seqfil=='N15NOE_lek_pfg_sel_enh_600'):
            zcoor='ncyc'


        elif(seqfil=='tnnoesy'):
            zcoor='mix'

        elif(seqfil=='tnnoesy_ajb'):
            zcoor='mix'



        elif(seqfil=='N15T2_lek_pfg_sel_enh_600'):
            zcoor='ncyc'
            pwn=float(vpar.GetParVarian('./procpar','n',('','pwn'))[0])
            delayFac=(32.0*pwn*1E-6 + 32.0*450.0e-6)

        elif(seqfil=='N15T1_lek_pfg_sel_enh_600'):
            zcoor='ncyc'
            pw_shpss=float(vpar.GetParVarian('./procpar','n',('','pw_shpss'))[0])
            delayFac = (pw_shpss*1E-6 + 2.0*2.5e-3);

        elif(seqfil=='NHT1_unenhanced_lek_600'):
            zcoor='ncyc'
            delayFac = (2.0*12.5e-3);





#        elif(seqfil=='CH3_1HT2s_600_lek'):
#            ConvertSpec_ProtonT2s(ni=ni)
#        elif(seqfil=='CH3_T1Z_T1ZZ_lek_600'):
#            ConvertSpec_IzIzz(ni=ni)
#        elif(seqfil=='hmqc_c13_600_methyl_lek'):
#            ConvertSpec_chmqc(ni=ni,ppmMin=ppmMin,ppmMax=ppmMax)
#        elif(seqfil=='hsqc_gd_sl_seduce_600'):
#            ConvertSpec_nhsqc(ni=ni,ppmMin=ppmMin,ppmMax=ppmMax)
#        elif(seqfil=='CT_hsqc_600'):
#            ConvertSpec_ct_chsqc(ni=ni,ppmMin=ppmMin,ppmMax=ppmMax)
        elif(seqfil=='N15_CPMG_Rex_NH_trosy_antitrosy_lek_600_v4'):
            zcoor='ncyc'
        else:
            print 'Cannot find ',seqfil,' in varian processing library'
            sys.exit(100)
    if(type=='omega'):
        
        parfile=vpar.GetOmegaParFile()
        test=vpar.GetParOmega(parfile,'n',('','seq_source',))[0].split('/')
        seqfil=test[len(test)-1]

        print 'Sequence file:',seqfil

 #       if(seqfil=='ajb_hmqc_c13_methyl.s'):
 #           ConvertSpec_chmqc_omeg(parfile,ni=ni)
 #       elif(seqfil=='ajb_HtoC_CH3_exchange_600_ILV.s'):
 #           ConvertSpec_HtoC_omeg(parfile,ni=ni)
 #       elif(seqfil=='ajb_HtoC_CH3_exchange_750_ILV.s'):
 #           ConvertSpec_HtoC_omeg(parfile,ni=ni)
 #       elif(seqfil=='ajb_HtoC_CH3_exchange_950_ILV.s'):
 #           ConvertSpec_HtoC_omeg(parfile,ni=ni)
 #       else:
 #           print 'Cannot find ',seqfil,' in varian processing library'
 #           sys.exit(100)

    if(type=='bruk'):
        seqfil=vpar.GetParBruk('acqu','n',('','PULPROG',))[0].split('<')[1].split('>')[0]
        if(seqfil=='r2_disp_exp_3D.mk'):
             zcoor='ncyc'

    outy=open(fudafile,'w')
    # FOR FITTING 40oC SINGLE QUANTUM aB-crystallin DATA
    # Read peaklist and spectrum info
    # 
    outy.write('PEAKLIST=%s\n' % peaklist)
    outy.write('SPECFILE=%s\n' % outfile)
    outy.write('NOISE=3263.0\n')
    outy.write('ZCOOR=%s\n' % (zcoor))
    outy.write('DELAYFACTOR=%f\n'% (delayFac))
    outy.write('BASELINE=N\n')
    outy.write('VERBOSELEVEL=5\n')
    outy.write('PRINTDATA=Y\n')
    outy.write('LM=(MAXFEV=50;TOL=1e-3)\n')
    outy.write('#DISCARD_SLICES=(1)\n')
    outy.write('#BASELINE=Y\n')
    outy.write('#\n')
    outy.write('#Specify the default values. All values are in ppm:\n')
    outy.write('#\n')
    outy.write('DEF_LINEWIDTH_F1=%f\n' % (radiusC*fac))
    outy.write('DEF_LINEWIDTH_F2=%f\n' % (radiusH*fac))
    outy.write('DEF_RADIUS_F1=%f\n' % radiusC )
    outy.write('DEF_RADIUS_F2=%f\n' % radiusH )
    outy.write('SHAPE=GLORE\n')
    outy.write('ISOTOPESHIFT=N\n')
    outy.write('#\n')
    outy.write('##\n')
    outy.write('######\n')
    outy.write('#DUMPPARAMETERS=Y\n')

    from baldwinStd import readfile
    array=readfile(peaklist)
    spec=[]
    for i in range(len(array)):
        try:
            name=array[i][0]
            if(float(array[i][2])>float(array[i][1])):
                cppm=float(array[i][2])
                hppm=float(array[i][1])
            else:
                cppm=float(array[i][1])
                hppm=float(array[i][2])
            spec.append((name,cppm,hppm))
        except:
            pass
        
    groups=FindOverlapGroups(spec,radiusC,radiusH)

    for i in range(len(groups)):
        outy.write('OVERLAP_PEAKS=(')
        for j in range(len(groups[i])):
            if(j!=0):
                outy.write(';')
            outy.write('%s' % (groups[i][j]))
        outy.write(')\n')
    outy.close()


def GetNeighbours(i,spec,radiusC,radiusH,grp,grpno):
    maxyC=spec[i][1]+radiusC
    minyC=spec[i][1]-radiusC
    maxyH=spec[i][2]+radiusH
    minyH=spec[i][2]-radiusH
    #print minyC,maxyC,minyH,maxyH
    for j in range(len(spec)):
        if(i!=j):
            if(spec[j][1]<maxyC and spec[j][1]>minyC):
                if(spec[j][2]<maxyH and spec[j][2]>minyH):
                    val=AddIfNew(spec[j][0],grp)
                    if(val==1):
                        grpno.append(j)
    return grp,grpno




def FindOverlapGroups(spec,radiusC,radiusH):
    overlap=[]
    for i in range(len(spec)):
        grp=[]
        grpno=[]
        go=0
        grp.append(spec[i][0])
        grpno.append(i)
        grp,grpno=GetNeighbours(i,spec,radiusC,radiusH,grp,grpno)
        while(go==0):
            start=len(grp)
            for j in range(len(grp)):
                grp,grpno=GetNeighbours(grpno[j],spec,radiusC,radiusH,grp,grpno)
            if(len(grp)==start):#finished
                go=1
        if(len(grp)>1):#if the group is complete
            AddIfNew(sorted(grp),overlap)
    return overlap





