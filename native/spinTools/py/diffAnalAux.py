#!/usr/bin/python

#######################################################################
# Class for diffusion analysis
#


import os,sys,baldwinStd,string,stats,numpy,baldwinStd
from scipy.optimize import leastsq,curve_fit
from scipy.special import sici
from numpy import arccos,sin,cos,log,exp,histogram2d,linspace
from scipy.special import betainc
import vpar,nmrglue as ng,math



class DiffyAnal():
    def __init__(self,pth,exp='n'):
        self.T=300;             #thermodynamics temperature (K)
        self.k=1.381E-23        #boltzmann's constant (Js-1)
        self.nu=0.001           #viscosity of water (cp)
        self.Gmax=60.           #G cm-1
        self.gamma=2.675222E4;  #gamma rad s-1 G-1
        self.conv=0.002         #conversion of DAC to Gcm-1
        
        #ppms to extract raw data for baselining and display
        self.ppmMin=-1   
        self.ppmMax=11. #12
        
        #ppm range to extract for global analysis
        self.analMin=0.82
        self.analMax=0.95


        self.ONLY2='n'
        #ppm range from which to estimate noise
        self.noiseMin=-1.
        self.noiseMax=-0.5
        self.boot=100  # number of runs when doing the boostraps
        self.datapath='./raw'  #location of the raw_{d20 in ms}.txt files eg raw_300.txt

        self.monexp_fit='y'  #run monexp fits on all shifts
        self.biexp_fit='y'   #run biexp fits on all shifts

        self.IntMax='*'
        self.rodmodel='n'

        self.base_lin='n'
        self.base_const='n'
        self.baseMin=-1.0
        self.baseMax=0.0


        self.pth=pth
        
        self.MonLim=3 #min number of points for a monomer fit


        if(exp=='n'):
            self.seqfil=vpar.GetParVarian(self.pth+'/procpar','n',('','seqfil'))[0].split('\"')[1]

            if(self.seqfil=="water_sLED_fm_v2_600" or self.seqfil=="water_sLED_fm_v2_600_19F" ):
                self.d20=vpar.GetParVarian(self.pth+'/procpar','n',('','BigT'))  #get d2o
                self.gradsRaw=vpar.GetParVarian(self.pth+'/procpar','n',('','gzlvl1')) #get gradients
                self.delta=float(vpar.GetParVarian(self.pth+'/procpar','n',('','gt1'))[0])*1000. # in ms

            elif(self.seqfil=="CT_N_hsqc_LED_lek_600_v2"):
                self.d20=vpar.GetParVarian(self.pth+'/procpar','n',('','Big_delta'))
                self.delta=float(vpar.GetParVarian(self.pth+'/procpar','n',('','gt2'))[0])*1000
                self.gradsRaw=vpar.GetParVarian(self.pth+'/procpar','n',('','gzlvl2'))
            elif(self.seqfil=="hmqc_c13_600_methyl_diffusion_lek"):
                self.d20=vpar.GetParVarian(self.pth+'/procpar','n',('','bigT'))
                self.delta=2*float(vpar.GetParVarian(self.pth+'/procpar','n',('','gt5'))[0])*1000
                self.gradsRaw=vpar.GetParVarian(self.pth+'/procpar','n',('','gzlvl5'))                   
        else:
            #self.d20=vpar.GetParBruk(self.pth+'/procpar','n',('','Big_delta'))            
            self.gradsRaw=[]
            for i in range(len(exp)):
                self.d20=vpar.GetBrukVal(exp[i]+'/acqus','D20')[0] #in s

                self.delta=float(vpar.GetBrukVal(exp[i]+'/acqus','P30')[0])/1000. #in ms
                print exp[i],'bigDelta: ', self.d20,' littleDelta: ',self.delta


                self.gradsRaw.append(float(vpar.GetBrukVal(exp[i]+'/acqus','GPZ6')[0])) #in ms
            print self.gradsRaw



    def FuncCalc(self,A0,grads,b,Deff):
        return A0*numpy.exp(-grads*b*Deff)


    def FitFuncRot(self,x,b):
        a=x[0]
        
        grads=b[0]
        data=b[1]
        d20=b[2]
        
        delta=b[3]
        Ycalc= self.Roddiff(a,grads,d20,delta)
        return Ycalc-data

    def FitFuncRotMix(self,x,b):
        a=numpy.fabs(x[0])
        Deff=numpy.fabs(x[1])
        mix=numpy.fabs(x[2])
        grads=b[0]
        data=b[1]
        d20=b[2]
        delta=b[3]
        bfac=(delta*self.gamma*self.Gmax)**2.*(d20-delta) #delays are already in seconds
        Ycalc= (1-mix)*self.Roddiff(a,grads,d20,delta)+self.FuncCalc(mix,grads,bfac,Deff)
        return Ycalc-data


    def FitFuncDiff(self,x,b):
        Deff=numpy.fabs(x[0])
        A0=numpy.fabs(x[1])
        grads=b[0]
        data=b[1]
        b=b[2]
        Ycalc=self.FuncCalc(A0,grads,b,Deff)
        return Ycalc-data

    def FitFuncDiffBi(self,x,b):
        Deff1=numpy.fabs(x[0])
        A01=numpy.fabs(x[2])
        Deff2=numpy.fabs(x[1])
        A02=numpy.fabs(x[3])
        grads=b[0]
        data=b[1]
        b=b[2]
        Ycalc=self.FuncCalc(A01,grads,b,Deff1)+self.FuncCalc(A02,grads,b,Deff2)
        return Ycalc-data

    #Estimate diffusion coefficient for initial calc
    def EstDeff(self,grads,inty,bfac):
        maxVal=grads[0]
        imax=0
        for i in range(len(grads)):
            if(grads[i]>maxVal):
                maxVal=grads[i]
                imax=i
        minVal=grads[0]
        imin=0
        for i in range(len(grads)):
            if(grads[i]<minVal):
                minVal=grads[i]
                imin=i
        #y1=A0*numpy.exp(-grads1*b*Deff)
        #y2=A0*numpy.exp(-grads2*b*Deff)
        #log(y2/y1)/(grads2/grads1)/b=Deff
        Deff=numpy.log(inty[imax]/inty[imin])/(-maxVal+minVal)/bfac
        return Deff


    def Roddiff(self,a,G,d20,delta):

        b0=1E-9;  #rod radius (m)
        
        N=(d20-delta/3)*self.Gmax*self.Gmax*delta*delta*self.gamma*self.gamma; 
        alpha = 1E2*delta*self.Gmax*self.gamma*pow(G,0.5); #10^2 puts into G into m-1 from cm-1
        frod = a/(log(a/b0)-0.3); #friction for trans diff
        d=(self.k*self.T)/(3*numpy.pi*self.nu*frod);   #trans diff
    
        #test=0  #for freely rotating
        #test=1  #for correlated motion (needs rotational diffusion)
        #test=2  #for zero rotational contribution
        test=0
        
        r4 = exp(-1*N*G*d*1E4);  #10^4 conv diff coeff into cm-2 needed for N
        
        if(test==0):
            r2=a;                                 
            zed=alpha*a;          
            sinseries1=sici(zed)[0];
            r1=(cos(zed)-1.+zed*sinseries1)/(0.5*zed*alpha);   #freely rotating
            yfree=r4*r1/r2
            return yfree

        if(test==1):
            r2=a;                                 
            frrod= (a*a*a)/(3*log(a/(2*b0))); #friction for rot diff
            dr=(k*T/(numpy.pi*nu*frrod));  #rotational diffusion coeff
            rr=(1-exp(-2*dr*d20));   #restricted rotation factor
            rr5=(rr**0.5);           
            zed2=alpha*a*rr5;
            sinseries2=sici(zed2)[0];
            r5=(cos(zed2)-1+zed2*sinseries2)/(0.5*zed2*alpha*rr5); #correlated
            ycorr=r4*r5/r2
            return ycorr

        if(test==2):
            ystat=r4
            return ystat



    def rawEdit(self,din):

        dinNew=[]
        shifts=[]
        for i in range(len(din)):
            if(i==0):
                grads=[]
                for j in range(len(din[i])-1):        
                    try:
                        test=float(din[i][j+1])
                        grads.append(float(din[i][j+1]))
                    except:
                        pass
                print 'Detected ',len(grads),'gradient strengths'
                #line=[]
                #line.append('#')
                #for i in range(len(grads)):
                #    line.append(grads[i])
                #dinNew.append(line)
            else:
                if(len(din[i])>len(grads)):
                    line=[]
                    shifts.append(float(din[i][1]))
                    for j in range(len(grads)):        
                        line.append(float(din[i][j+2]))
                    dinNew.append(line)
        return numpy.array(grads),numpy.array(shifts),numpy.array(dinNew)


    #extract shift range and apply linear baseline correction
    def shiftExtract(self,din,shifts):

        shiftsNew=[]
        dinNew=[]
        for i in range(len(shifts)):
            if(shifts[i]>self.ppmMin and shifts[i]<self.ppmMax):
                shiftsNew.append(shifts[i])
                dinNew.append(din[i])
        spec=[]
        dinNewNew=numpy.empty_like(dinNew)

        #subtract a straightline from the start to the end point

        if(self.base_lin=='y'):
            for j in range(len(dinNew[0])):
                x1=shiftsNew[0]
                y1=dinNew[0][j]
                x2=shiftsNew[len(shiftsNew)-1]
                y2=dinNew[len(dinNew)-1][j]
        
                #y2=m*x2+c
                #y1=m*x1+c
            
                m=(y2-y1)/(x2-x1)
                c=y2-m*x2
        
                for i in range(len(dinNew)):
                    dinNewNew[i][j]=dinNew[i][j]-(m*shiftsNew[i]+c)
        elif(self.base_const=='y'):
            for j in range(len(dinNew[0])):
                shiftsNew=numpy.array(shiftsNew)
                dinNew=numpy.array(dinNew)
                
                mask=(shiftsNew>self.baseMin)*(shiftsNew<self.baseMax) 
                print dinNew[mask],self.baseMin,self.baseMax
                c=numpy.average(dinNew[mask])
                
                for i in range(len(dinNew)):
                    dinNewNew[i][j]=dinNew[i][j]-(c)
        else:
            dinNewNew=dinNew


        return numpy.array(shiftsNew),numpy.array(dinNewNew)


    def MakeData(self,infile,bruk='n'):
        if(bruk=='n'):
            d20=self.d20
            gradsRaw=self.gradsRaw

            grads=[]
            for i in range(len(gradsRaw)):
                grads.append( (float(gradsRaw[i])*self.conv/self.Gmax*100.))  #gradients in per cent
            print 'Gradients (%): ',grads
    
            dic,data = ng.pipe.read(self.pth+'/'+infile)
            Size=data.shape
            uc0 = ng.pipe.make_uc(dic,data,dim=0)
            uc1 = ng.pipe.make_uc(dic,data,dim=1)
        
            uc0max=uc0.ppm(0)
            uc0min=uc0.ppm(Size[0]-1)
            uc1max=uc1.ppm(0)
            uc1min=uc1.ppm(Size[1]-1)
        
            print "Spectrum dimensions (pts): ",Size   #print the spectral dimensions
            print "dimension 0 limits (ppm): ", uc0min, uc0max  #carbon 
            print "dimension 1 limits (ppm): ", uc1min, uc1max  #direct 
        
            baldwinStd.PathExists((self.datapath,)) #create path if not there already
            #output as a single 'raw' file
            for k in range(len(d20)):
                outy=open(self.datapath+'/raw_'+str(int(float(d20[k])*1000))+'.txt','w')
                outy.write('No.\tChemical shift \t')
                for i in range(len(grads)):
                    outy.write('%s\t' % grads[i])
                outy.write('\n')
                for i in range(Size[1]):
                    outy.write('%i %e\t' % (i,uc1.ppm(i)))
                    for j in range(len(grads)):
                        outy.write('%e\t' % data[j+k*len(grads),i])
                    outy.write('\n')
                outy.close()

        else:
            d20=self.d20,
            gradsRaw=self.gradsRaw

            grads=[]
            for i in range(len(gradsRaw)):
                grads.append( (float(gradsRaw[i]) ))  #gradients in per cent
            print 'Gradients (%): ',grads
    
            data=[]
            for i in range(len(infile)):
                print infile[i]
                dicA,dataA = ng.pipe.read(infile[i]+'/test.ft2')
                Size=dataA.shape
                uc0 = ng.pipe.make_uc(dicA,dataA,dim=0)
                #uc1 = ng.pipe.make_uc(dic,data,dim=1)
                data.append(dataA)
                uc0max=uc0.ppm(0)
                uc0min=uc0.ppm(Size[0]-1)
                #uc1max=uc1.ppm(0)
                #uc1min=uc1.ppm(Size[1]-1)
        
                print "Spectrum dimensions (pts): ",Size   #print the spectral dimensions
                print "dimension 0 limits (ppm): ", uc0min, uc0max  #carbon 
                #print "dimension 1 limits (ppm): ", uc1min, uc1max  #direct 
        
            baldwinStd.PathExists((self.datapath,)) #create path if not there already

            #output as a single 'raw' file
            for k in range(len(d20)):
                outy=open(self.datapath+'/raw_'+str(int(float(d20[k])*1000))+'.txt','w')
                outy.write('No.\tChemical shift \t')
                for i in range(len(grads)):
                    outy.write('%s\t' % grads[i])
                outy.write('\n')
                for i in range(Size[0]):
                    outy.write('%i %e\t' % (i,uc0.ppm(i)))
                    for j in range(len(grads)):
                        outy.write('%e\t' % data[j][i])
                    outy.write('\n')
                outy.close()



        self.GetData()


    def GetData(self): #get BigT values of available files
        filey=os.listdir(self.datapath)
        d20=[]
        for i in range(len(filey)):#search for valid data files
            if(filey[i][0:3]=='raw' and filey[i][len(filey[i])-4:]=='.txt'):
                d20val=float(string.split(filey[i],'_')[1].split('.')[0])
                d20.append(d20val)
        self.d20=sorted(d20)


    def DiffAnal(self,grads,shifts,din,d20val):

        noiseVals=[]
        for i in range(len(shifts)):
            if(shifts[i]>self.noiseMin and shifts[i]<self.noiseMax):
                for j in range(len(din[i])):
                    noiseVals.append(din[i][j])
        if(len(noiseVals)==0):
            print 'No noise! aborting.'
            sys.exit(100)
        #print noiseVals
        noise=max(noiseVals)
        print 'Noise:',noise

        
        outy=open(self.datapath+'/data_'+str(int(d20val))+'.txt','w')
        outy.write('#\t')
        for j in range(len(grads)):
            outy.write('%f\t' % ((grads[j])))    
        outy.write('\n')

        Shuf=[]
        Diff=[]
        Aval=[]
        
        Shuf2=[]
        Diff1=[]
        Diff2=[]
        Aval1=[]
        Aval2=[]

        for i in range(len(din)):
            ex=[]
            ey=[]
            for j in range(len(din[i])):
                if(din[i][j]>noise):
                    ex.append((grads[j]/100)**2.) #factor of Gmax^2
                    ey.append(din[i][j])
            ex=numpy.array(ex)
            ey=numpy.array(ey)
                
            #        grads.append(((float(din[0][j+1])/100.)**2.))

            if(len(ey)>self.MonLim):
                if(self.monexp_fit=='y'):
                    bfac=(self.delta*1E-3*self.gamma*self.Gmax)**2.*(d20val-self.delta)*1E-3 #set conversion factor

                    Deff=self.EstDeff(ex,ey,bfac)
                    x0=leastsq(self.FitFuncDiff,[Deff,max(ey)],args=[ex,ey,bfac])
          
                    Shuf.append(shifts[i])
                    Diff.append(math.fabs(x0[0][0]))
                    Aval.append(x0[0][1])

                    if(shifts[i]>self.analMin and shifts[i]<self.analMax):
                        outy.write('%f\t' % (shifts[i]))
                        for j in range(len(grads)):
                            outy.write('%f\t' % (din[i][j]/x0[0][1]))
                        outy.write('\n')

                if(len(ey)>6):

                    if(self.biexp_fit=='y'):
                        x0=leastsq(self.FitFuncDiffBi,[Deff/5,Deff*5,max(ey)/2.,max(ey)/2.,],args=[ex,ey,bfac])

                        Shuf2.append(shifts[i])
                        Diff1.append(math.fabs(x0[0][0]))
                        Aval1.append(math.fabs(x0[0][2]))
                        Diff2.append(math.fabs(x0[0][1]))
                        Aval2.append(math.fabs(x0[0][3]))
                else:
                    Shuf2.append(shifts[i])
                    Aval1.append(Aval[len(Aval)-1])
                    Diff1.append(Diff[len(Diff)-1])
                    Aval2.append(0.0)
                    Diff2.append(1E-12)
                    

                #if(shifts[i]>analMin and shifts[i]<analMax):
                #    outy.write('%f\t' % (shifts[i]))
                #    for j in range(len(grads)):
                #        outy.write('%f\t' % (din[i][j]/x0[0][1]))
                #    outy.write('\n')
            
            else:
                Shuf2.append(shifts[i])
                Aval1.append(0.0)
                Diff1.append(1E-12)
                Aval2.append(0.0)
                Diff2.append(1E-12)           

 
        outy.close()

        #print the diffusion data to a file
        if(self.monexp_fit=='y'):
            outy=open(self.datapath+'/diff_'+str(int(d20val))+'.txt','w')
            for j in range(len(Shuf)):
                outy.write('%f\t%e\t%e\n' % (Shuf[j],Diff[j],Aval[j]))
            outy.write('\n\n')
            if(self.biexp_fit=='y'):
                for j in range(len(Shuf2)):
                    outy.write('%f\t%e\t%e\t%e\t%e\n' % (Shuf2[j],Diff1[j],Diff2[j],Aval1[j],Aval2[j]))
            outy.close()

            #make histograms
            print 'Making histograms'

            if(numpy.amin(Shuf)>0):
                binX=linspace(numpy.amin(Shuf)*0.95,numpy.amax(Shuf)*1.05,self.binsPPM)
            else:
                binX=linspace(numpy.amin(Shuf)*1.05,numpy.amax(Shuf)*0.95,self.binsPPM)

            Diff=numpy.fabs(Diff)
            Diff1=numpy.array(Diff1)
            Diff2=numpy.array(Diff2)


            DiffMax=numpy.amax(Diff)*1.1
            DiffMax1=numpy.amax(Diff1)*1.1
            DiffMax2=numpy.amax(Diff2)*1.1

            if(DiffMax>DiffMax1):
                DiffMax=DiffMax
            else:
                DiffMax=DiffMax1

            if(DiffMax1>DiffMax2):
                DiffMax=DiffMax1
            else:
                DiffMax=DiffMax2

            self.HistMax=DiffMax
            DiffMin=self.HistMin #set minimum value from histogram values


            binY= DiffMin*10.**(numpy.log10(DiffMax/DiffMin)*(numpy.arange(self.binsDiff)/(self.binsDiff-1.)))


            H,xedges,yedges=histogram2d(Shuf,Diff,bins=(binX,binY))
            H1,xedges,yedges=histogram2d(Shuf2,Diff1,bins=(binX,binY))
            H2,xedges,yedges=histogram2d(Shuf2,Diff2,bins=(binX,binY))

            outy=open(self.datapath+'/histy_'+str(int(d20val))+'.out','w')               
            for i in range(len(binX)-1):
                for j in range(len(binY)-1):
                    outy.write('%e\t%e\t%f\n' % ((binX[i]+binX[i+1])/2.,(binY[j]+binY[j+1])/2.,H[i,j]))
                outy.write('\n')
            outy.write('\n\n')
            for i in range(len(binX)-1):
                for j in range(len(binY)-1):
                    outy.write('%e\t%e\t%f\n' % ((binX[i]+binX[i+1])/2.,(binY[j]+binY[j+1])/2.,H1[i,j]+H2[i,j]))
                outy.write('\n')            
            outy.close()

            outy=open(self.datapath+'/histyProj_'+str(int(d20val))+'.out','w')               
            for j in range(len(binY)-1):
                outy.write('%e\t%f\n' % ((binY[j]+binY[j+1])/2.,numpy.sum(H[:,j])))
            outy.write('\n\n')
            for j in range(len(binY)-1):
                outy.write('%e\t%f\n' % ((binY[j]+binY[j+1])/2.,numpy.sum(H1[:,j]+H2[:,j])))
            outy.close()




            

            thresh=self.sepLine
            
            mk1=Diff1>thresh
            mk2=Diff1<DiffMax
            mk3=Diff2>thresh
            mk4=Diff2<DiffMax
            Asmall=Aval1*mk1*mk2+Aval2*mk3*mk4


            mk1=Diff1>DiffMin
            mk2=Diff1<thresh
            mk3=Diff2>DiffMin
            mk4=Diff2<thresh
            Abig=Aval1*mk1*mk2+Aval2*mk3*mk4
            
            outy=open(self.datapath+'/decon_'+str(int(d20val))+'.out','w')
            for i in range(len(Shuf2)):
                outy.write('%e\t%e\t%e\n' % (Shuf2[i],Asmall[i],Abig[i]))          
                #if(Shuf2[i+1]-Shuf2[i]>shifts[1]-shifts[0]):
                #    outy.write('%e\t%e\t%e\n' % (Shuf2[i],Asmall[i],Abig[i]))          
            outy.close()
        return 


    def ReadProc(self,infile):
    
        array=baldwinStd.readfile(infile)
        
        dinNew=[]
        shifts=[]
        for i in range(len(array)):
            if(i==0):
                grads=[]
                for j in range(len(array[i])-1):        
                    try:
                        test=float(array[i][j+1])
                        grads.append(float(array[i][j+1]))
                    except:
                        pass
                print 'Detected ',len(grads),'gradient strengths'
                #line=[]
                #line.append('#')
                #for i in range(len(grads)):
                #    line.append(grads[i])
                #dinNew.append(line)
            else:
                if(len(array[i])>len(grads)):
                    line=[]
                    shifts.append(float(array[i][0]))
                    for j in range(len(grads)):        
                        line.append(float(array[i][j+1]))
                    dinNew.append(line)
        return numpy.array(grads),numpy.array(shifts),numpy.array(dinNew)


    def ProcData(self):
        for i in range(len(self.d20)):
            print 'Reading in raw data file: '+self.datapath+'/raw_'+str(int(self.d20[i]))+'.txt'
            din=baldwinStd.readfile(self.datapath+'/raw_'+str(int(self.d20[i]))+'.txt')
            
            grads,shifts,din=self.rawEdit(din) #edit dataformat

            shifts,din=self.shiftExtract(din,shifts) #extract specified region and baseline

            #print out processed file
            outy=open(self.datapath+'/proc_'+str(int(self.d20[i]))+'.txt','w')
            outy.write('#\t')
            for j in range(len(grads)):
                outy.write('%f\t' % (grads[j]))
            outy.write('\n')
            for j in range(len(din)):
                outy.write('%f\t' % (shifts[j]))
                for k in range(len(din[j])):
                    outy.write('%f\t' % (din[j][k]))
                outy.write('\n')
            outy.close()

    def AnalProc(self):
        for i in range(len(self.d20)):
            print 'Reading in processed data file: '+self.datapath+'/proc_'+str(int(self.d20[i]))+'.txt'
            grads,shifts,din=self.ReadProc(self.datapath+'/proc_'+str(int(self.d20[i]))+'.txt')    
            #perform diffusion analysis of each spectrum
            #extract the anal region to a data*.txt file for subsequent analysis
            self.DiffAnal(grads,shifts,din,self.d20[i]) #only take data above noise


            #make gnuplot script to print the spectra, with the value exrapolated from the fit
            baldwinStd.PathExists(('gnu',))
            gnu=open('gnu/spec.gp','w')
            gnu.write('set term post eps enh color solid\n')
            gnu.write('set output \''+self.datapath+'/raw_'+str(int(self.d20[i]))+'.txt.eps\'\n')
            gnu.write('set xlabel \'1H (ppm)\'\n')
            gnu.write('set title \''+str(int(self.d20[i]))+'ms raw data\'\n')
            gnu.write('plot ')
            for j in range(len(grads)):
                gnu.write('\''+self.datapath+'/proc_%s.txt\' u 1:%i ti \'%s\' w li lt %i' %(str(int(self.d20[i])),j+2,grads[j],j+2))
                #if(i!=len(grads)-1):
                gnu.write(',')
            gnu.write('\''+self.datapath+'/diff_%s.txt\' i 0 u 1:3 ti \'%s\' w li lt %i' %(str(int(self.d20[i])),'extrapolated (mon exp)',1))
            gnu.write('\n')    
            
            gnu.write('set output \''+self.datapath+'/raw_'+str(int(self.d20[i]))+'.1.txt.eps\'\n')
            gnu.write('set yrange [*:%s]\n' % (str(self.IntMax)))
            gnu.write('set title \''+str(int(self.d20[i]))+'ms (zoomed intensity)\'\n')
            gnu.write('plot ')
            for j in range(len(grads)):
                gnu.write('\''+self.datapath+'/proc_%s.txt\' u 1:%i ti \'%s\' w li lt %i' %(str(int(self.d20[i])),j+2,grads[j],j+2))
                #if(i!=len(grads)-1):
                gnu.write(',')
            gnu.write('\''+self.datapath+'/diff_%s.txt\' i 0 u 1:3 ti \'%s\' w li lt %i' %(str(int(self.d20[i])),'extrapolated (mon exp)',1))
            gnu.write('\n')    

            gnu.write('set output \''+self.datapath+'/diff_'+str(int(self.d20[i]))+'.eps\'\n')
            gnu.write('set logscale y\n')
            gnu.write('plot \''+self.datapath+'/diff_%s.txt\' i 0 u 1:2 ti \'monoexp\'w li,\'\' i 1 u 1:2 ti \'biexp1\',\'\' i 1 u 1:3 ti \'biexp2\'' %(str(int(self.d20[i]))))
            gnu.write('\n')    


            gnu.close()
            os.system('gnuplot gnu/spec.gp') #run gnuplot script


        #overlay all diffusion data from all diffusion times
        gnu=open('gnu/diff.gp','w')
        gnu.write('set term post eps enh color solid\n')
        gnu.write('set output \''+self.datapath+'/diff.eps\'\n')
        gnu.write('set logscale y\n')
        gnu.write('set xlabel \'1H (ppm)\'\n')
        gnu.write('set ylabel \'D_{eff} (cm^2s^{-1})\'\n')
        gnu.write('set format y "%.0t*10^%T"\n')
        gnu.write('set yrange[%s:%s]\n' % (str(self.HistMin),str(self.HistMax)))
        gnu.write('set title \'Monoexponential diffusion fits\'\n')
        gnu.write('plot ')
        for i in range(len(self.d20)):
            gnu.write('\''+self.datapath+'/diff_%s.txt\' i 0 u 1:2 ti \'%s\' lt %i' %(str(int(self.d20[i])),str(int(self.d20[i])),i+1))
            if(i!=len(self.d20)-1):
                gnu.write(',')
        gnu.write('\n')    

        if(self.ONLY2=='n'):
            gnu.write('set output \''+self.datapath+'/diffBi.eps\'\n')
            gnu.write('set logscale y\n')
            gnu.write('set xlabel \'1H (ppm)\'\n')
            gnu.write('set ylabel \'D_{eff} (cm^2s^{-1})\'\n')
            gnu.write('set title \'Biexponential diffusion fits\'\n')
            gnu.write('set yrange[%s:%s]\n' % (str(self.HistMin),str(self.HistMax)))
            gnu.write('plot ')
            for i in range(len(self.d20)):
                gnu.write('\''+self.datapath+'/diff_%s.txt\' i 1 u 1:2 ti \'%s\' lt %i,\'\' i 1 u 1:3 noti lt %i' %(str(int(self.d20[i])),str(int(self.d20[i])),i+1,i+1))
                if(i!=len(self.d20)-1):
                    gnu.write(',')
            gnu.write('\n')    
    


        gnu.write('set output \''+self.datapath+'/diffDecon.eps\'\n')
        gnu.write('set title \'Reconstructed spectra above and below %.2e cm^2s^{-1}\'\n' % (self.sepLine))
        gnu.write('set key\n')
        gnu.write('set yrange[*:*]\n')
        gnu.write('unset logscale\n')
        gnu.write('set xlabel \'1H (ppm)\'\n')
        gnu.write('set ylabel \'I\'\n')
        gnu.write('plot ')
        gnu.write('\''+self.datapath+'/decon_%s.out\' u 1:2 ti \'big species %s\' w li,\\\n' %(str(int(self.d20[0])),str(int(self.d20[0]))))
        gnu.write('\''+self.datapath+'/decon_%s.out\' u 1:3 ti \'small species %s\' w li\n' %(str(int(self.d20[0])),str(int(self.d20[0]))))

        gnu.write('set output \''+self.datapath+'/diffDecon1.eps\'\n')
        gnu.write('set title \'Reconstructed spectra above and below %.2e cm^2s^{-1}\'\n' % (self.sepLine))
        gnu.write('set key\n')
        gnu.write('set yrange[*:%s]\n' % (str(self.IntMax)))
        gnu.write('unset logscale\n')
        gnu.write('set xlabel \'1H (ppm)\'\n')
        gnu.write('set ylabel \'I\'\n')
        gnu.write('plot ')
        gnu.write('\''+self.datapath+'/decon_%s.out\' u 1:2 ti \'big species %s\' w li,\\\n' %(str(int(self.d20[0])),str(int(self.d20[0]))))
        gnu.write('\''+self.datapath+'/decon_%s.out\' u 1:3 ti \'small species %s\' w li\n' %(str(int(self.d20[0])),str(int(self.d20[0]))))



        gnu.write('set output \''+self.datapath+'/diffHistProj1.eps\'\n')
        gnu.write('set logscale x\n')
        gnu.write('set key\n')
        gnu.write('set yrange[*:*]\n')
        gnu.write('set title \'Diffusion probability histograms (ppm projections)\'\n')
        gnu.write('set ylabel\'Probabilty\'\n')
        gnu.write('set xlabel \'D_{eff} (cm^2s^{-1})\'\n')
        gnu.write('set format x "%.0t*10^%T"\n')
        gnu.write('plot ')
        gnu.write('\''+self.datapath+'/histyProj_%s.out\' i 0 u 1:2 ti \'monoExp %s\' w li,\\\n' %(str(int(self.d20[0])),str(int(self.d20[0]))))
        gnu.write('\''+self.datapath+'/histyProj_%s.out\' i 1 u 1:2 ti \'biExp %s\' w li\n' %(str(int(self.d20[0])),str(int(self.d20[0]))))
        gnu.write('\n')    



        gnu.write('set output \''+self.datapath+'/diffHist1.eps\'\n')
        gnu.write('set logscale y\n')
        gnu.write('unset logscale x\n')
        gnu.write('set format x "%f"\n')
        gnu.write('set cntrparam levels incr 0,5,100\n' )
        gnu.write('unset key\n')
        gnu.write('set xlabel \'1H (ppm)\'\n')
        gnu.write('set ylabel \'D_{eff} (cm^2s^{-1})\'\n')
        gnu.write('set format y "%.0t*10^%T"\n')
        gnu.write('set view map\n')
        gnu.write('unset surface\n')
        gnu.write('set contour\n')
        gnu.write('splot ')
        gnu.write('\''+self.datapath+'/histy_%s.out\' i 0 u 1:2:3 ti \'%s\' w li lc palette\n' %(str(int(self.d20[0])),str(int(self.d20[0]))))
        gnu.write('\n')    


        gnu.write('set output \''+self.datapath+'/diffHist2.eps\'\n')
        gnu.write('unset key\n')
        gnu.write('set logscale y\n')
        gnu.write('set xlabel \'1H (ppm)\'\n')
        gnu.write('set ylabel \'D_{eff} (cm^2s^{-1})\'\n')
        gnu.write('set format y "%.0t*10^%T"\n')
        gnu.write('set view map\n')
        gnu.write('unset surface\n')
        gnu.write('set contour\n')
        gnu.write('splot ')
        gnu.write('\''+self.datapath+'/histy_%s.out\' i 1 u 1:2:3 ti \'%s\' w li lc palette\n' %(str(int(self.d20[0])),str(int(self.d20[0]))))
        gnu.write('\n')    
    

        gnu.close()
        os.system('gnuplot gnu/diff.gp')

        if(self.ONLY2=='n'):
            os.system('arraygraph.py 4 8 0 0 0 0 `ls '+self.datapath+'/raw*.eps` '+self.datapath+'/diff.eps '+self.datapath+'/diffBi.eps '+self.datapath+'/diffHist1.eps '+self.datapath+'/diffHist2.eps '+self.datapath+'/diffDecon1.eps '+self.datapath+'/diffDecon.eps '+self.datapath+'/diffHistProj1.eps ')
        else:
            os.system('arraygraph.py 4 8 0 0 0 0 `ls '+self.datapath+'/raw*.eps` '+self.datapath+'/diff.eps '+self.datapath+'/diffHist1.eps '+self.datapath+'/diffHist2.eps '+self.datapath+'/diffDecon1.eps '+self.datapath+'/diffDecon.eps '+self.datapath+'/diffHistProj1.eps ')


        os.system('mv summary.pdf rawData.pdf')
    




    def AnalData(self):
        
        baldwinStd.PathExists(('out',))

        #array for calculated gradient strengths
        Gcalc=numpy.linspace(self.Gcmin,self.Gcmax,self.Gcrid)
        Gcalc=(Gcalc/100.) **2.


        outy=open('out/test.out','w')
        outy2=open('out/calc.out','w')
        outy3=open('out/param.out','w')

        sele=[]
        #sele=2,5,9  #pick which traces, numbered from 0
        #sele=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14
        
        dataFull=[]
        for i in range(len(self.d20)):#for each d20
            print 'reading in '+self.datapath+'/data_'+str(int(self.d20[i]))+'.txt'
            
            d20val=self.d20[i]
            din=baldwinStd.readfile(self.datapath+'/data_'+str(int(self.d20[i]))+'.txt')
            
            grads=[]#get gradient squared in units of ratio
            for j in range(len(din[0])-1):
                grads.append(((float(din[0][j+1])/100.)**2.))
            grads=numpy.array(grads)


            #get intensities and standard deviations of decay curves
            inty=[]
            err=[]
            for j in range(len(grads)):
                vals=[]
                if(len(sele)!=0):
                    for k in range(len(sele)):
                        vals.append(float(din[sele[k]+1][j+1]))
                else:
                    for k in range(len(din)-1):
                        vals.append(float(din[k+1][j+1]))
                        
                inty.append((stats.mean(vals) ))
                err.append((stats.stdev(vals) ))

                
            inty=numpy.array(inty) #integrate
            err=numpy.array(err)  #integrate

            dataFull.append((grads,inty))

            #run LM best fit to averaged data
            bfac=(self.delta*1E-3*self.gamma*self.Gmax)**2.*(d20val-self.delta)*1E-3 #set conversion factor
            Deff=self.EstDeff(grads,inty,bfac)
            x0=leastsq(self.FitFuncDiff,[Deff,max(inty)],args=[grads,inty,bfac])
            
            print 'Bestfit: ',x0[0][0]
            #run a boot strap to get numerical errors
            Dboot=[]
            for j in range(self.boot):
                gradsNew=[]
                intyNew=[]
                for k in range(len(grads)):#select grads points with replacement
                    test=numpy.random.randint(0,len(grads))
                    gradsNew.append(grads[test])
                    intyNew.append(inty[test])
                gradsNew=numpy.array(gradsNew)
                intyNew=numpy.array(intyNew)
                Deff=self.EstDeff(gradsNew,intyNew,bfac)
                x0new=leastsq(self.FitFuncDiff,[Deff,max(intyNew)],args=[gradsNew,intyNew,bfac])
                Dboot.append(x0new[0][0])

            print 'Bootstrap mean:  ',stats.mean(Dboot)
            print 'Bootstrap stdev: ',stats.stdev(Dboot)

            try:
                x0_2=leastsq(self.FitFuncDiffBi,[x0[0][0]/5,x0[0][0]*5,max(inty)/2.,max(inty)/2.,],args=[grads,inty,bfac])
        
                #print params
                outy3.write('%f\t%e\t%f\t%e\t%e\t%e\t%e\t%e\t%e\n' % (d20val,x0[0][0],x0[0][1],stats.mean(Dboot),stats.stdev(Dboot),math.fabs(x0_2[0][0]),math.fabs(x0_2[0][1]),math.fabs(x0_2[0][2]),math.fabs(x0_2[0][3])))
            except:
                pass

            #print raw data to text file
            for j in range(len(grads)):
                outy.write('%f\t%s\t%s\n' % ((grads[j],inty[j],err[j])))
            outy.write('\n\n')

            #print fitted data
            calc=self.FuncCalc(x0[0][1],Gcalc,bfac,x0[0][0])

            for j in range(int(self.Gcrid)):
                outy2.write('%f\t%f\n' % ((Gcalc[j],calc[j])))
            outy2.write('\n\n')
            try:
                calc2=self.FuncCalc(math.fabs(x0_2[0][2]),Gcalc,bfac,math.fabs(x0_2[0][0]))+self.FuncCalc(math.fabs(x0_2[0][3]),Gcalc,bfac,math.fabs(x0_2[0][1]))
                for j in range(int(self.Gcrid)):
                    outy2.write('%f\t%f\n' % ((Gcalc[j],calc2[j])))
            except:
                pass
            outy2.write('\n\n')

        outy.close()
        outy2.close()
        outy3.close()
        self.dataFull=dataFull





    def ProcUntangle(self,din):
        dinNew=[]
        grads=[]
        shifts=[]
        
        for i in range(len(din)):
            if(i==0):
                for j in range(len(din[0])-1):
                    grads.append(float(din[0][j+1]))
            else:
                shifts.append(float(din[i][0]))
                line=[]
                for j in range(len(din[i])-1):
                    line.append(float(din[i][j+1]))
                if(len(line) != len(grads)):
                    print 'error!'
                    sys.exit(100)
                dinNew.append(line)
        return numpy.array(shifts),numpy.array(grads),numpy.array(dinNew)


    def FitSpecInt(self,x,b):
        dataFull=b[0]
        d20Full=b[1]
        gradsFull=b[2]
        shiftsFull=b[3]
        bfacFull=b[4]
        delta=b[5]
        d20N=b[6]
        shiftsN=b[7]
        gradsN=b[8]
        mode=b[9]

        if(mode==0):
            ycalc=self.CalcInty1(x,d20N,shiftsN,gradsN,d20Full,delta,gradsFull,bfacFull)
        elif(mode==1):
            ycalc=self.CalcInty2(x,d20N,shiftsN,gradsN,d20Full,delta,gradsFull,bfacFull)
        elif(mode==2):
            ycalc=self.CalcInty3(x,d20N,shiftsN,gradsN,d20Full,delta,gradsFull,bfacFull)

        #    return (ycalc-dataFull)*gradsFull*10000.
        return (ycalc-dataFull)



    def CalcInty(self,x,d20N,shiftsN,gradsN,rod,d20Full,delta,mix,gradsFull,bfacFull,Deff):
        A0=numpy.empty_like(gradsFull)
        for i in range(d20N):#unpack intensities
            for k in range(shiftsN):
                for j in range(gradsN):
                    A0[j+k*gradsN+i*shiftsN*gradsN]=x[k+i*shiftsN]
        return A0*((1-mix)*self.Roddiff(rod,gradsFull,d20Full,delta*1E-3)+self.FuncCalc(mix,gradsFull,bfacFull,Deff))

    def CalcInty1(self,x,d20N,shiftsN,gradsN,d20Full,delta,gradsFull,bfacFull):
        #get first spectrum
        spec1=numpy.empty(shiftsN)
        for i in range(shiftsN):
            spec1[i]=x[i]
        spec1=spec1/numpy.max(spec1)

        #get multipliers
        d20mult=numpy.empty(d20N)
        for i in range(d20N):
            d20mult[i]=x[shiftsN+i]
                  
        #get Deff
        Deff=x[shiftsN+d20N]

        #get intensity multipliers
        A0=numpy.empty_like(gradsFull)
        for i in range(d20N):#unpack intensities
            for k in range(shiftsN):
                for j in range(gradsN):
                    A0[j+k*gradsN+i*shiftsN*gradsN]=d20mult[i]*spec1[k]
                    
        return self.FuncCalc(A0,gradsFull,bfacFull,Deff)


    def CalcInty2(self,x,d20N,shiftsN,gradsN,d20Full,delta,gradsFull,bfacFull):
        #get first spectrum
        spec1=numpy.empty(shiftsN)
        for i in range(shiftsN):
            spec1[i]=x[i]
        spec1=numpy.abs(spec1)
        spec1=spec1/numpy.max(spec1)

        #get second spectrum
        spec2=numpy.empty(shiftsN)
        for i in range(shiftsN):
            spec2[i]=x[i+shiftsN]
        spec2=numpy.abs(spec2)
        spec2=spec2/numpy.max(spec2)

        
        d20mult=numpy.empty(d20N)
        for i in range(d20N):
            d20mult[i]=x[shiftsN*2+i]

        Deff=x[shiftsN*2+d20N]
        Deff2=x[shiftsN*2+d20N+1]
        mix=x[shiftsN*2+d20N+2]

        #get intensity multipliers
        A0rod=numpy.empty_like(gradsFull)
        for i in range(d20N):#unpack intensities
            for k in range(shiftsN):
                for j in range(gradsN):
                    A0rod[j+k*gradsN+i*shiftsN*gradsN]=d20mult[i]*spec2[k]

        #get intensity multipliers
        A0free=numpy.empty_like(gradsFull)
        for i in range(d20N):#unpack intensities
            for k in range(shiftsN):
                for j in range(gradsN):
                    A0free[j+k*gradsN+i*shiftsN*gradsN]=d20mult[i]*spec1[k]

        return A0rod*(self.FuncCalc((1.-mix),gradsFull,bfacFull,Deff))+A0free*(self.FuncCalc(mix,gradsFull,bfacFull,Deff2))


    def CalcInty3(self,x,d20N,shiftsN,gradsN,d20Full,delta,gradsFull,bfacFull):
        #get first spectrum
        spec1=numpy.empty(shiftsN)
        for i in range(shiftsN):
            spec1[i]=x[i]
        spec1=spec1/numpy.max(spec1)

        #get second spectrum
        spec2=numpy.empty(shiftsN)
        for i in range(shiftsN):
            spec2[i]=x[i+shiftsN]
        spec2=spec2/numpy.max(spec2)

        d20mult=numpy.empty(d20N)
        for i in range(d20N):
            d20mult[i]=x[shiftsN*2+i]
            
        rod=x[shiftsN*2+d20N]
        mix=x[shiftsN*2+d20N+1]
        Deff=x[shiftsN*2+d20N+2]

        #get intensity multipliers
        A0rod=numpy.empty_like(gradsFull)
        for i in range(d20N):#unpack intensities
            for k in range(shiftsN):
                for j in range(gradsN):
                    A0rod[j+k*gradsN+i*shiftsN*gradsN]=d20mult[i]*spec2[k]

        #get intensity multipliers
        A0free=numpy.empty_like(gradsFull)
        for i in range(d20N):#unpack intensities
            for k in range(shiftsN):
                for j in range(gradsN):
                    A0free[j+k*gradsN+i*shiftsN*gradsN]=d20mult[i]*spec1[k]
                    
        return A0rod*((1-mix)*self.Roddiff(rod,gradsFull,d20Full,delta*1E-3))+A0free*(self.FuncCalc(mix,gradsFull,bfacFull,Deff))


    def SmoothData(self,smooth,shifts,grads,din):
        
        specFull=[]
        cnt=0
        shiAve=[]
        intAve=[]
        
        for j in range(len(grads)):
            specNew=[]
            shiftsNew=[]
            for k in range(len(shifts)):
                intAve.append(din[k][j])
                shiAve.append(shifts[k])
            
                cnt+=1
                if(cnt==smooth):
                    shiftsNew.append(numpy.average(shiAve))
                    specNew.append(numpy.average(intAve))
                    intAve=[]
                    shiAve=[]
                    cnt=0
            cnt=0 #reset the count
            intAve=[]
            shiAve=[]
            specNew=numpy.array(specNew)
            specFull.append(specNew)

        #transpose specFull
        dinNew=[]
        for k in range(len(shiftsNew)):
            line=[]
            for j in range(len(grads)):
                line.append(specFull[j][k])
            dinNew.append(line)
        dinNew=numpy.array(dinNew)
        return numpy.array(shiftsNew),numpy.array(dinNew)/numpy.max(dinNew)

    ##############################################################################
    #Functions for global fits

    def FitRodMon(self,smooth,d20,shifts,grads,dataFull,d20Full,shiftsFull,gradsFull,specMax,mix,rod,Deff):
        #SETUP ROD+monomer fit

        #take the first spectrum as the initial template for both species
        #normalise to 1
        specTemp1=numpy.empty_like(shifts)
        specTemp2=numpy.empty_like(shifts)
        for k in range(len(shifts)):
            specTemp1[k]=dataFull[0+k*len(grads)+0*len(shifts)*len(grads)]
            specTemp2[k]=dataFull[0+k*len(grads)+0*len(shifts)*len(grads)]

        specTemp1=specTemp1/numpy.max(specTemp1)
        specTemp2=specTemp2/numpy.max(specTemp2)

        #construct array containing minmisation parameters
        initVal=[]
        for i in range(len(specTemp1)):
            initVal.append(specTemp1[i] )
        for i in range(len(specTemp2)):
            initVal.append(specTemp2[i] )
        for i in range(len(specMax)):
            initVal.append(specMax[i])
        initVal.append(rod)
        initVal.append(mix)
        initVal.append(Deff)

        initVal=numpy.array(initVal)
        
        bfacFull=(delta*1E-3*gamma*Gmax)**2.*(d20Full-delta*1E-3) 
        bfacFull=numpy.array(bfacFull)

        min_flg='y'
        if(min_flg=='y'):
            mode=2 #set to rod mode
            print 'Performing global fit to rod+monomer model (slow)'
            print 'NOTE: - increasing smooth to speed up.'
            print 'Current set to ',smooth,' with: '
            print 'Points per spectrum:    ',len(shifts)
            print 'Gradient points per d2o:',len(grads)
            print 'D20 values:             ',len(d20)
            print 'Total datapoints:       ',len(d20)*len(grads)*len(shifts)
            x0=leastsq(self.FitSpecInt,initVal,args=[dataFull,d20Full,gradsFull,shiftsFull,bfacFull,self.delta,len(d20),len(shifts),len(grads),mode],full_output=1)
            print 'Minimisation complete'
        else:
            x0=initVal,

        print 'd20 scale vals:'
        for i in range(len(d20)):
            print d20[i],x0[0][len(shifts)*2+i]

        print 'Rod before:',initVal[len(shifts)*2+len(d20)]
        rod=x0[0][len(shifts)*2+len(d20)]
        print 'Rod after:',rod
        print 'mix before:',initVal[len(shifts)*2+len(d20)+1]
        mix=x0[0][len(shifts)*2+len(d20)+1]
        print 'mix after:',mix
        print 'Deff before:',initVal[len(shifts)*2+len(d20)+2]
        Deff=x0[0][len(shifts)*2+len(d20)+2]
        print 'Deff after:',Deff

        if(min_flg=='y'):
            roderr=x0[1][len(shifts)*2+len(d20),len(shifts)*2+len(d20)]**0.5
            mixerr=x0[1][len(shifts)*2+len(d20)+1,len(shifts)*2+len(d20)+1]**0.5
            Defferr=x0[1][len(shifts)*2+len(d20)+2,len(shifts)*2+len(d20)+2]**0.5
        
        os.system('rm glob/*.eps')

        ycalc=self.CalcInty3(initVal,len(d20),len(shifts),len(grads),d20Full,delta,gradsFull,bfacFull)
        chi2=((ycalc-dataFull)**2.).sum()
        print 'chi2:',chi2
        print 'red chi2:',chi2/(len(dataFull)-len(x0[0]))

        #print fitted data file
        outy=open('glob/sim.out','w')
        for j in range(len(grads)):
            for k in range(len(shifts)):
                outy.write('%e\t' % shiftsFull[j+k*len(grads)])
                for i in range(len(d20)):
                    outy.write('%e\t' % (ycalc[j+k*len(grads)+i*len(grads)*len(shifts)]))
                outy.write('\n')
            outy.write('\n\n')
        outy.close()
    

        #update the specs
        for i in range(len(shifts)):
            specTemp1[i]=x0[0][i]
            specTemp2[i]=x0[0][i+len(shifts)]
        specTemp1=specTemp1/numpy.max(specTemp1)
        specTemp2=specTemp2/numpy.max(specTemp2)

        outy=open('glob/calc.out','w')
        for i in range(len(shifts)):
            outy.write('%e\t%e\t%e\t%e\t%e\n' % (shifts[i],mix*(specTemp1[i]),(1-mix)*(specTemp2[i]),specTemp1[i]/numpy.max(specTemp1),specTemp2[i]/numpy.max(specTemp2)))
        outy.close()

        gnu=open('gnu/fit.gp','w')
        gnu.write('set term post eps enh color solid\n')
        for i in range(len(d20)):
            gnu.write('set output \'glob/fit.'+str(int(d20[i]))+'.eps\'\n')
            gnu.write('set xlabel \'1H (ppm)\'\n')
            gnu.write('set title \'fitted spectrum: '+str(int(d20[i]))+'ms\'\n')
            gnu.write('plot ')
            gnu.write('\'glob/data.out\' u 1:%i noti w li lt 1,\'glob/sim.out\' u 1:%i noti w li lt 2' %(i+2,i+2))
            if(i!=len(d20)-1):
                gnu.write(',')
            gnu.write('\n')    

        gnu.write('set output \'glob/decon.eps\'\n')
        if(min_flg=='y'):
            gnu.write('set label sprintf(\"rod: %s.2f +/- %s.2f\",%f,%f) at graph 0.6,0.9\n' % ('%','%',1-mix,mixerr))
            gnu.write('set label sprintf(\"free: %s.2f +/- %s.2f\",%f,%f) at graph 0.6,0.85\n' % ('%','%',mix,mixerr))
            gnu.write('set label sprintf(\"Monomer Deff: %s.2e +/- %s.2e\",%e,%e) at graph 0.6,0.8\n' % ('%','%',Deff,Defferr))
            gnu.write('set label sprintf(\"Rod length: %s.2e +/- %s.2e\",%e,%e) at graph 0.6,0.75\n' % ('%','%',rod,roderr))
        else:
            gnu.write('set label sprintf(\"rod: %s.2f\",%f) at graph 0.7,0.9\n' % ('%',1-mix))
            gnu.write('set label sprintf(\"free: %s.2f\",%f) at graph 0.7,0.85\n' % ('%',mix))
            gnu.write('set label sprintf(\"Monomer Deff: %s.2e\",%e) at graph 0.7,0.8\n' % ('%',Deff))
            gnu.write('set label sprintf(\"Rod length: %s.2e\",%e) at graph 0.7,0.75\n' % ('%',rod))
        gnu.write('set label sprintf(\"chi2: %s.2e\",%e) at graph 0.7,0.7\n' % ('%',chi2))
        gnu.write('set label sprintf(\"red chi2: %s.2e\",%e) at graph 0.7,0.65\n' % ('%',chi2/(len(dataFull)-len(x0[0]))))
            
        gnu.write('set title \'Deconvolved spectra at accurate relative intensities\'\n')
        gnu.write('plot \'glob/calc.out\' u 1:2 ti \'free species\' w li lt 1,\'\' u 1:3 ti \'rod\' w li lt 2')
        gnu.write('\n')    
        gnu.write('set output \'glob/decon2.eps\'\n')
        gnu.write('set title \'Deconvolved spectra at normalised intensities\'\n')
        gnu.write('plot \'glob/calc.out\' u 1:4 ti \'free species\' w li lt 1,\'\' u 1:5 ti \'rod\' w li lt 2')
        
        gnu.close()
        os.system('gnuplot gnu/fit.gp') #run gnuplot script
        
        os.system('arraygraph.py 4 8 0 0 0 0 `ls glob/*.eps`')
        os.system('mv summary.pdf fit.RodMon.pdf')
        
        return chi2,len(dataFull)-len(x0[0])



    def FitTwoMon(self,smooth,d20,shifts,grads,dataFull,d20Full,shiftsFull,gradsFull,specMax,Deff,mix):

        #SETUP monomer fit
        specTemp1=numpy.empty_like(shifts)
        for k in range(len(shifts)):
            specTemp1[k]=dataFull[0+k*len(grads)+0*len(shifts)*len(grads)]
        specTemp2=numpy.empty_like(shifts)
        for k in range(len(shifts)):
            specTemp2[k]=dataFull[0+k*len(grads)+0*len(shifts)*len(grads)]

        initVal=[]
        for i in range(len(specTemp1)):
            initVal.append(specTemp1[i] )
        for i in range(len(specTemp2)):
            initVal.append(specTemp2[i] )
        for i in range(len(specMax)):
            initVal.append(specMax[i])
        initVal.append(Deff)
        initVal.append(Deff/100.)
        initVal.append(mix)
        initVal=numpy.array(initVal)
        
        bfacFull=(self.delta*1E-3*self.gamma*self.Gmax)**2.*(d20Full-self.delta*1E-3) 
        bfacFull=numpy.array(bfacFull)
        min_flg='y'
        if(min_flg=='y'):
            mode=1 #set to two monomer  fit
            print 'Performing global fit to two monomer model'
            print 'NOTE: - increasing smooth to speed up.'
            print 'Current set to ',smooth,' with: '
            print 'Points per spectrum:    ',len(shifts)
            print 'Gradient points per d2o:',len(grads)
            print 'D20 values:             ',len(d20)
            print 'Total datapoints:       ',len(d20)*len(grads)*len(shifts)
            x0=leastsq(self.FitSpecInt,initVal,args=[dataFull,d20Full,gradsFull,shiftsFull,bfacFull,self.delta,len(d20),len(shifts),len(grads),mode],full_output=1)
            print 'Minimisation complete'
        else:
            x0=initVal,

        print 'd20 scale vals:'
        for i in range(len(d20)):
            print d20[i],x0[0][len(shifts)*2+i]

        print 'Deff before:',initVal[len(shifts)*2+len(d20)]
        Deff=x0[0][len(shifts)*2+len(d20)]
        print 'Deff after:',Deff
        print 'Deff2 before:',initVal[len(shifts)*2+len(d20)+1]
        Deff2=x0[0][len(shifts)*2+len(d20)+1]
        print 'Deff2 after:',Deff2
        print 'mix before:',initVal[len(shifts)*2+len(d20)+2]
        mix=x0[0][len(shifts)*2+len(d20)+2]
        print 'mix after:',mix

        if(min_flg=='y'):
            try:
                Defferr=x0[1][len(shifts)*2+len(d20),len(shifts)*2+len(d20)]**0.5
                Deff2err=x0[1][len(shifts)*2+len(d20)+1,len(shifts)*2+len(d20)+1]**0.5
                mixerr=x0[1][len(shifts)*2+len(d20)+2,len(shifts)*2+len(d20)+2]**0.5
            except:
                Defferr=1.
                Deff2err=1.
                mixerr=1.

        ycalc=self.CalcInty2(x0[0],len(d20),len(shifts),len(grads),d20Full,self.delta,gradsFull,bfacFull)

        chi2=((ycalc-dataFull)**2.).sum()
        print 'chi2:',chi2
        print 'red chi2:',chi2/(len(dataFull)-len(x0[0]))
        
        os.system('rm glob/*.eps')

        #print fitted data file
        outy=open('glob/sim.out','w')
        for j in range(len(grads)):
            for k in range(len(shifts)):
                outy.write('%e\t' % shiftsFull[j+k*len(grads)])
                for i in range(len(d20)):
                    outy.write('%e\t' % (ycalc[j+k*len(grads)+i*len(grads)*len(shifts)]))
                outy.write('\n')
            outy.write('\n\n')
        outy.close()
    

        #update the specs
        for i in range(len(shifts)):
            specTemp1[i]=x0[0][i]
        specTemp1=numpy.abs(specTemp1)
        specTemp1=specTemp1/numpy.max(specTemp1)
        for i in range(len(shifts)):
            specTemp2[i]=x0[0][i+len(shifts)]
        specTemp2=numpy.abs(specTemp2)
        specTemp2=specTemp2/numpy.max(specTemp2)

        outy=open('glob/calc.out','w')
        for i in range(len(shifts)):
            outy.write('%e\t%e\t%e\t%e\t%e\n' % (shifts[i],mix*(specTemp1[i]),(1-mix)*specTemp2[i],specTemp1[i]/numpy.max(specTemp1),specTemp2[i]/numpy.max(specTemp2)))
        outy.close()

        gnu=open('gnu/fit.gp','w')
        gnu.write('set term post eps enh color solid\n')
        for i in range(len(d20)):
            gnu.write('set output \'glob/fit.'+str(int(d20[i]))+'.eps\'\n')
            gnu.write('set xlabel \'1H (ppm)\'\n')
            gnu.write('set title \'fitted spectrum: '+str(int(d20[i]))+'ms\'\n')
            gnu.write('plot ')
            gnu.write('\'glob/data.out\' u 1:%i noti w li lt 1,\'glob/sim.out\' u 1:%i noti w li lt 2' %(i+2,i+2))
            if(i!=len(d20)-1):
                gnu.write(',')
            gnu.write('\n')    

        gnu.write('set output \'glob/decon.eps\'\n')
        if(min_flg=='y'):

            gnu.write('set label sprintf(\"Monomer Deff: %s.2e +/- %s.2e\",%e,%e) at graph 0.6,0.9\n' % ('%','%',Deff,Defferr))
            gnu.write('set label sprintf(\"Monomer Deff2: %s.2e +/- %s.2e\",%e,%e) at graph 0.6,0.85\n' % ('%','%',Deff2,Deff2err))
            gnu.write('set label sprintf(\"mix: %s.2e +/- %s.2e\",%e,%e) at graph 0.6,0.8\n' % ('%','%',mix,mixerr))
        
        else:
            gnu.write('set label sprintf(\"Monomer Deff: %s.2e\",%e) at graph 0.7,0.9\n' % ('%',Deff))
            gnu.write('set label sprintf(\"Monomer Deff2: %s.2e\",%e) at graph 0.7,0.85\n' % ('%',Deff2))
            gnu.write('set label sprintf(\"mix: %s.2e\",%e) at graph 0.7,0.8\n' % ('%',mix))
        gnu.write('set label sprintf(\"chi2: %s.2e\",%e) at graph 0.7,0.75\n' % ('%',chi2))
        gnu.write('set label sprintf(\"red chi2: %s.2e\",%e) at graph 0.7,0.7\n' % ('%',chi2/(len(dataFull)-len(x0[0]))))

        gnu.write('set title \'Deconvolved spectra at accurate relative intensities\'\n')
        gnu.write('plot \'glob/calc.out\' u 1:2 ti \'free species1\' w li lt 1,\'\' u 1:3 ti \'free species2\' w li lt 2')
        gnu.write('\n') 
        gnu.write('set output \'glob/decon2.eps\'\n')   
        gnu.write('set title \'Deconvolved spectra at normalised relative intensities\'\n')
        gnu.write('plot \'glob/calc.out\' u 1:4 ti \'free species1\' w li lt 1,\'\' u 1:5 ti \'free species2\' w li lt 2')

        gnu.close()
        os.system('gnuplot gnu/fit.gp') #run gnuplot script
        
        os.system('arraygraph.py 4 8 0 0 0 0 `ls glob/*.eps`')
        os.system('mv summary.pdf fit.MonMon.pdf')
        

        return chi2,len(dataFull)-len(x0[0])





    def FitMonOnly(self,smooth,d20,shifts,grads,dataFull,d20Full,shiftsFull,gradsFull,specMax,Deff):

        #SETUP monomer fit
        specTemp1=numpy.empty_like(shifts)
        for k in range(len(shifts)):
            specTemp1[k]=dataFull[0+k*len(grads)+0*len(shifts)*len(grads)]
        initVal=[]
        for i in range(len(specTemp1)):
            initVal.append(specTemp1[i] )
        for i in range(len(specMax)):
            initVal.append(specMax[i])
        initVal.append(Deff)
        initVal=numpy.array(initVal)

        bfacFull=(self.delta*1E-3*self.gamma*self.Gmax)**2.*(d20Full-self.delta*1E-3) 
        bfacFull=numpy.array(bfacFull)
        min_flg='y'
        if(min_flg=='y'):
            mode=0 #set to monomer only fit
            print 'Performing global fit to monomer model'
            print 'NOTE: - increasing smooth to speed up.'
            print 'Current set to ',smooth,' with: '
            print 'Points per spectrum:    ',len(shifts)
            print 'Gradient points per d2o:',len(grads)
            print 'D20 values:             ',len(d20)
            print 'Total datapoints:       ',len(d20)*len(grads)*len(shifts)
            x0=leastsq(self.FitSpecInt,initVal,args=[dataFull,d20Full,gradsFull,shiftsFull,bfacFull,self.delta,len(d20),len(shifts),len(grads),mode],full_output=1)
            print 'Minimisation complete'
        else:
            x0=initVal,

        print 'd20 scale vals:'
        for i in range(len(d20)):
            print d20[i],x0[0][len(shifts)+i]
            
        print 'Deff before:',Deff
        Deff=x0[0][len(shifts)+len(d20)]
        print 'Deff after:',Deff
        if(min_flg=='y'):
            Defferr=x0[1][len(shifts)+len(d20),len(shifts)+len(d20)]**0.5
    
        ycalc=self.CalcInty1(x0[0],len(d20),len(shifts),len(grads),d20Full,self.delta,gradsFull,bfacFull)

        chi2=((ycalc-dataFull)**2.).sum()
        print 'chi2:',chi2
        print 'red chi2:',chi2/(len(dataFull)-len(x0[0]))

        os.system('rm glob/*.eps')

        #print fitted data file
        outy=open('glob/sim.out','w')
        for j in range(len(grads)):
            for k in range(len(shifts)):
                outy.write('%e\t' % shiftsFull[j+k*len(grads)])
                for i in range(len(d20)):
                    outy.write('%e\t' % (ycalc[j+k*len(grads)+i*len(grads)*len(shifts)]))
                outy.write('\n')
            outy.write('\n\n')
        outy.close()
    

        #update the specs
        for i in range(len(shifts)):
            specTemp1[i]=x0[0][i]
        specTemp1=specTemp1/numpy.max(specTemp1)

        outy=open('glob/calc.out','w')
        for i in range(len(shifts)):
            outy.write('%e\t%e\t%e\n' % (shifts[i],(specTemp1[i]),specTemp1[i]/numpy.max(specTemp1)))
        outy.close()

        gnu=open('gnu/fit.gp','w')
        gnu.write('set term post eps enh color solid\n')
        for i in range(len(d20)):
            gnu.write('set output \'glob/fit.'+str(int(d20[i]))+'.eps\'\n')
            gnu.write('set xlabel \'1H (ppm)\'\n')
            gnu.write('set title \'fitted spectrum: '+str(int(d20[i]))+'ms\'\n')
            gnu.write('plot ')
            gnu.write('\'glob/data.out\' u 1:%i noti w li lt 1,\'glob/sim.out\' u 1:%i noti w li lt 2' %(i+2,i+2))
            if(i!=len(d20)-1):
                gnu.write(',')
            gnu.write('\n')    

        gnu.write('set output \'glob/decon.eps\'\n')
        if(min_flg=='y'):
            gnu.write('set label sprintf(\"Monomer Deff: %s.2e +/- %s.2e\",%e,%e) at graph 0.6,0.9\n' % ('%','%',Deff,Defferr))
        else:
            gnu.write('set label sprintf(\"Monomer Deff: %s.2e\",%e) at graph 0.7,0.9\n' % ('%',Deff))
        gnu.write('set label sprintf(\"chi2: %s.2e\",%e) at graph 0.7,0.85\n' % ('%',chi2))
        gnu.write('set label sprintf(\"red chi2: %s.2e\",%e) at graph 0.7,0.8\n' % ('%',chi2/(len(dataFull)-len(x0[0]))))

        gnu.write('set title \'Deconvolved spectra at accurate relative intensities\'\n')
        gnu.write('plot \'glob/calc.out\' u 1:2 ti \'free species\' w li lt 1')
        gnu.write('\n')    

        gnu.close()
        os.system('gnuplot gnu/fit.gp') #run gnuplot script
        
        os.system('arraygraph.py 4 8 0 0 0 0 `ls glob/*.eps`')
        os.system('mv summary.pdf fit.MonOnly.pdf')
        
        return chi2,len(dataFull)-len(x0[0])






    def SpecDecon(self):
        
        baldwinStd.PathExists(('glob',))
        
        d20Full=[]
        gradsFull=[]
        dataFull=[]
        shiftsFull=[]
        specMax=[]
        for i in range(len(self.d20)):
            print 'Reading in processed data file: '+self.datapath+'/proc_'+str(int(self.d20[i]))+'.txt'
            din=baldwinStd.readfile(self.datapath+'/proc_'+str(int(self.d20[i]))+'.txt')
        
            shifts,grads,din=self.ProcUntangle(din)

            smooth=self.smooth
            shifts,din=self.SmoothData(smooth,shifts,grads,din)


            tmpVal=[] #store all intensities for all gradients/shifts
            for k in range(len(shifts)):
                tmpVal.append(max(din[k]))
                for j in range(len(grads)):
                    shiftsFull.append((shifts[k]))
                    gradsFull.append(((grads[j])/100.)**2.)
                    d20Full.append((self.d20[i])*1E-3)
                    dataFull.append(din[k][j])
            specMax.append(max(tmpVal))#the most intense for all gradients/shifts


        dataFull=numpy.array(dataFull)
        d20Full=numpy.array(d20Full)
        gradsFull=numpy.array(gradsFull)
        shiftsFull=numpy.array(shiftsFull)

        #print data files
        outy=open('glob/data.out','w')
        for j in range(len(grads)):
            for k in range(len(shifts)):
                outy.write('%e\t' % shiftsFull[j+k*len(grads)])
                for i in range(len(self.d20)):
                    outy.write('%e\t' % (dataFull[j+k*len(grads)+i*len(grads)*len(shifts)]))
                outy.write('\n')
            outy.write('\n\n')
        outy.close()

        chi2Mon,dofMon      =self.FitMonOnly(smooth,self.d20,shifts,grads,dataFull,d20Full,shiftsFull,gradsFull,specMax,self.Deff)
        if(self.ONLY2=='n'):
            chi2TwoMon,dofTwoMon=self.FitTwoMon(smooth,self.d20,shifts,grads,dataFull,d20Full,shiftsFull,gradsFull,specMax,self.Deff,self.mix)
        if(self.rodmodel=='y'):
            chi2RodMon,dofRodMon=self.FitRodMon(smooth,self.d20,shifts,grads,dataFull,d20Full,shiftsFull,gradsFull,specMax,self.mix,self.rod,self.Deff)

        if(self.ONLY2=='n'):
            Plvl1=self.GetPlvl(chi2Mon,chi2TwoMon,dofMon,dofTwoMon)
            sys.stdout.write('Plvl mon-twomon:    %.2e' % Plvl1)
        if(self.rodmodel=='y'):
            Plvl2=self.GetPlvl(chi2TwoMon,chi2RodMon,dofTwoMon,dofRodMon)
            Plvl3=self.GetPlvl(chi2Mon,chi2RodMon,dofMon,dofRodMon)
            sys.stdout.write('Plvl twomon-rodmon: %.2e' % Plvl2)
            sys.stdout.write('Plvl mon-rodmon:    %.2e' % Plvl3)


    def GetPlvl(self,chi2Simple,chi2Complex,dofSimple,dofComplex):
        v1=dofSimple*1.
        v2=dofComplex*1.
        Z=v1-v2              #difference in dof
        F=(chi2Simple*1.-chi2Complex*1.)/((v1-v2)*chi2Complex*1./v2*1.) #the F statistic
        ex=v2/(v2+Z*F)
        plvl=betainc(v2/2,Z/2,ex) #integrate to get Plvl
        return plvl



    def GlobData(self):
        #array for calculated gradient strengths
        Gcalc=numpy.linspace(self.Gcmin,self.Gcmax,self.Gcrid)
        Gcalc=(Gcalc/100.) **2.

        print 'Assembling data for global fitting'
        #concatenate all data into one giant data array
        gFull=[]
        intFull=[]
        d20Full=[]
        for i in range(len(self.dataFull)):
            for j in range(len(self.dataFull[i][0])):
                gFull.append(self.dataFull[i][0][j])
                intFull.append(self.dataFull[i][1][j])
                d20Full.append(self.d20[i]/1000.)  #d20values in ms
        gFull=numpy.array(gFull)
        d20Full=numpy.array(d20Full)
        intFull=numpy.array(intFull)

        print intFull.shape,d20Full.shape,gFull.shape
    
        #print FitFuncRot([1E-6],[gFull,intFull,d20Full,delta*1E-3])

        #fit to freely rotating rod model
        x0=leastsq(self.FitFuncRot,[1E-6],args=[gFull,intFull,d20Full,self.delta*1E-3])
        a=x0[0]
        print 'Best fit rod length:',a
        ycalc=self.Roddiff(a,gFull,d20Full,self.delta*1E-3)
        print 'Ave%Err:',100*numpy.sqrt(((ycalc-intFull)**2.).sum()/(len(ycalc)-1.))
        print 'redchi2:',((ycalc-intFull)**2.).sum()/(len(intFull)-len(x0[0]))


        #get numerical errors from fit through bootstrap
        Dboot=[]
        for j in range(self.boot):
            gradsNew=[]
            intyNew=[]
            d20New=[]
            for k in range(len(gFull)):
                test=numpy.random.randint(0,len(gFull))
                gradsNew.append(gFull[test])
                intyNew.append(intFull[test])
                d20New.append(d20Full[test])
            gradsNew=numpy.array(gradsNew)
            intyNew=numpy.array(intyNew)
            d20New=numpy.array(d20New)
            #    Deff=EstDeff(gradsNew,intyNew,bfac)
            x0new=leastsq(self.FitFuncRot,[1E-6],args=[gradsNew,intyNew,d20New,self.delta*1E-3])
            Dboot.append(x0new[0])
        print 'Bootstrap mean:  ',stats.mean(Dboot)
        print 'Bootstrap stdev: ',stats.stdev(Dboot)

        try:
            x0_2=leastsq(self.FitFuncRotMix,[1E-6,1E-6,0.9],args=[gFull,intFull,d20Full,self.delta*1E-3])

            rod=numpy.fabs(x0_2[0][0])
            DeffMon=numpy.fabs(x0_2[0][1])
            mix=numpy.fabs(x0_2[0][2])

            print 'Mixture of free and fibril model:'
            print 'Length:           ',rod
            print 'Deff monomer:     ',DeffMon
            print 'Molefrac Monomer: ',mix
        except:
            pass

        #Get simulated data curves and Deff values
        gFull=[]
        d20Full=[]
        for i in range(len(self.dataFull)):
            for j in range(len(Gcalc)):
                gFull.append(Gcalc[j])
                d20Full.append(self.d20[i]/1000.)  #d20values in ms
        gFull=numpy.array(gFull)
        d20Full=numpy.array(d20Full)
        ycalc=self.Roddiff(a,gFull,d20Full,self.delta*1E-3)



        bfac=(self.delta*1E-3*self.gamma*self.Gmax)**2.*(d20Full-self.delta*1E-3) 

        try:
            ycalc2=(1-numpy.fabs(x0_2[0][2]))*self.Roddiff(numpy.fabs(x0_2[0][0]),gFull,d20Full,self.delta*1E-3)+self.FuncCalc(numpy.fabs(x0_2[0][2]),gFull,bfac,numpy.fabs(x0_2[0][1]))

            ycalc3=(1-numpy.fabs(x0_2[0][2]))*self.Roddiff(numpy.fabs(x0_2[0][0]),gFull,d20Full,self.delta*1E-3)
        except:
            pass



        outy=open('out/sim.out','w')
        outy2=open('out/simDeff.out','w')
        for j in range(len(self.d20)):
            gTmp=[]
            iTmp=[]
            lTmp=[]
            jTmp=[]
            for i in range(len(Gcalc)):
                outy.write('%f\t%e\t' % (Gcalc[i],ycalc[i+len(Gcalc)*j]))
                try:
                    outy.write('%e' % (ycalc2[i+len(Gcalc)*j]))
                except:
                    pass

                outy.write('\n')
                gTmp.append(Gcalc[i])
                iTmp.append(ycalc[i+len(Gcalc)*j])
                try:
                    lTmp.append(ycalc2[i+len(Gcalc)*j])
                    jTmp.append(ycalc3[i+len(Gcalc)*j])
                except:
                    pass
            gTmp=numpy.array(gTmp)
            iTmp=numpy.array(iTmp)
            #fit the calculated decay curves to model to get effective diffusion rates
            bfac=(self.delta*1E-3*self.gamma*self.Gmax)**2.*(self.d20[j]-self.delta)*1E-3 #set conversion factor
            Deff=self.EstDeff(gTmp,iTmp,bfac)
            x0=leastsq(self.FitFuncDiff,[Deff,max(iTmp)],args=[gTmp,iTmp,bfac])

            try:
                x0b=leastsq(self.FitFuncDiff,[Deff,max(iTmp)],args=[gTmp,lTmp,bfac])
                x0c=leastsq(self.FitFuncDiff,[Deff,max(iTmp)],args=[gTmp,jTmp,bfac])
            except:
                pass
            
            outy2.write('%f\t%e\t' % (self.d20[j],x0[0][0]))

            try:
                outy2.write('%e\t%e\t%e\t' % (x0b[0][0],x0c[0][0],numpy.fabs(x0_2[0][1])))
            except:
                pass
            outy2.write('\n')
            outy.write('\n\n')
        outy.close()
        outy2.close()

        try:
            self.rod=rod
            self.mix=mix
            self.Deff=DeffMon
        except:
            self.Deff=stats.mean(Dboot)
            self.mix=0.5

    def MakePlots(self):
        #make outputs
        gnu=open('gnu/gnu.gp','w')
        gnu.write('set term post eps enh color solid\n')
        gnu.write('set size square\n')


        gnu.write('set output \'out/plot0.eps\'\n')
        gnu.write('set title \'I/I0 for a  mon exponential fit\'\n')
        gnu.write('set xlabel \'1H(ppm)\'\n')
        gnu.write('plot ')
        jj=0
        for i in range(len(self.d20)):
            din=baldwinStd.readfile(self.datapath+'/data_'+str(int(self.d20[i]))+'.txt')
            grads=[]#get gradient squared in units of ratio
            for j in range(len(din[0])-1):
                grads.append(((float(din[0][j+1])/100.)**2.))
            grads=numpy.array(grads)
            for j in range(len(grads)):
                gnu.write('\'%s/data_%s.txt\' u 1:($%i) ti \'%s %s\' w li lt %i' %(self.datapath,str(int(self.d20[i])),2+j,str(self.d20[i]),str(grads[j]),j+1))
                jj+=1
                if(jj!=len(grads)*len(self.d20)):
                    gnu.write(',')
        gnu.write('\n')

                

        gnu.write('set output \'out/plot1.eps\'\n')
        gnu.write('set title \'I/I0 vs G/Gmax^2 with individual fits for Deff\'\n')
        gnu.write('set xlabel \'(G/Gmax)^2\'\n')
        gnu.write('plot ')
        for i in range(len(self.d20)):
            gnu.write('\'out/test.out\' i %i u 1:2:3 ti \'%s\' w err lt %i,\'out/calc.out\' i %i u 1:2 noti w li  lt %i,' %(i,self.d20[i],i+1,2*i,i+1))
            gnu.write('\'out/calc.out\' i %i u 1:2 noti w li  lt %i' %(2*i+1,i+1)) #put in biexp line
            
            if(i!=len(self.d20)-1):
                gnu.write(',')
        gnu.write('\n')
        gnu.write('set output \'out/plot1b.eps\'\n')
        gnu.write('set title \'I/I0 vs G/Gmax^2 with individual fits for Deff - log\'\n')
        gnu.write('set logscale y\n')
        gnu.write('set xlabel \'(G/Gmax)^2\'\n')
        gnu.write('plot ')
        for i in range(len(self.d20)):
            gnu.write('\'out/test.out\' i %i u 1:2:3 ti \'%s\' w err lt %i,\'out/calc.out\' i %i u 1:2 noti w li  lt %i,' %(i,self.d20[i],i+1,2*i,i+1))
            gnu.write('\'out/calc.out\' i %i u 1:2 noti w li  lt %i' %(2*i+1,i+1)) #put in biexp line
            
            if(i!=len(self.d20)-1):
                gnu.write(',')
        gnu.write('\n')

        gnu.write('set output \'out/plot2.eps\'\n')
        gnu.write('unset logscale y\n')
        gnu.write('set size square\n')
        gnu.write('set title \'Deff vs bigDelta\'\n')
        gnu.write('set xlabel \'bigD (ms)\'\n')
        gnu.write('set ylabel \'Deff (10^{-6} cm^2s^{-1}) \'\n')
        gnu.write('set xrange[*:%f]\n' % (max(self.d20)*1.05))
        gnu.write('set size square\n')
        gnu.write('unset key\n')
        
        #gnu.write('set format y "%.2t*10^{%T}"\n')
        gnu.write('plot ')
        #gnu.write('\'param.out\' u 1:($4*1E6):($5*1E6) w err,\'\' u 1:($2*1E6)')
        gnu.write('\'out/param.out\' u 1:($4*1E6):($5*1E6) w err,\'out/simDeff.out\' u 1:($2*1E6) w li,\'\' u 1:($3*1E6) w li,\'\' u 1:($4*1E6) w li,\'\' u 1:($5*1E6) w li,\\\n')
        gnu.write('\'out/param.out\' u 1:($6*1E6),\'\' u 1:($7*1E6)')
        
        gnu.write('\n')
        gnu.write('reset\n')
        gnu.write('set term post eps enh color solid\n')
        gnu.write('set title \'I/I0 vs G/Gmax^2 with global fit (rod only)\'\n')
        gnu.write('set size square\n')
        gnu.write('set output \'out/plot3.eps\'\n')
        gnu.write('set xlabel \'(G/Gmax)^2\'\n')
        gnu.write('plot ')
        for i in range(len(self.d20)):
            gnu.write('\'out/test.out\' i %i u 1:2:3 ti \'%s\' w err lt %i,\'out/sim.out\' i %i u 1:2 noti w li  lt %i' %(i,self.d20[i],i+1,i,i+1))
            if(i!=len(self.d20)-1):
                gnu.write(',')
        gnu.write('\n')
        gnu.write('set title \'I/I0 vs G/Gmax^2 with global fit (rod only) - log plot\'\n')
        gnu.write('set logscale y\n')
        gnu.write('set size square\n')
        gnu.write('set output \'out/plot3b.eps\'\n')
        gnu.write('set xlabel \'(G/Gmax)^2\'\n')
        gnu.write('plot ')
        for i in range(len(self.d20)):
            gnu.write('\'out/test.out\' i %i u 1:2:3 ti \'%s\' w err lt %i,\'out/sim.out\' i %i u 1:2 noti w li  lt %i' %(i,self.d20[i],i+1,i,i+1))
            if(i!=len(self.d20)-1):
                gnu.write(',')
        gnu.write('\n')
        gnu.write('set title \'I/I0 vs G/Gmax^2 with global fit (rod and free monomer)\'\n')
        gnu.write('unset logscale y\n')
        gnu.write('set size square\n')
        gnu.write('set output \'out/plot4.eps\'\n')
        gnu.write('set xlabel \'(G/Gmax)^2\'\n')
        gnu.write('plot ')
        for i in range(len(self.d20)):
            gnu.write('\'out/test.out\' i %i u 1:2:3 ti \'%s\' w err lt %i,\'out/sim.out\' i %i u 1:3 noti w li  lt %i' %(i,self.d20[i],i+1,i,i+1))
            if(i!=len(self.d20)-1):
                gnu.write(',')
        gnu.write('\n')
        gnu.write('set title \'I/I0 vs G/Gmax^2 with global fit (rod and free monomer) - log plot\'\n')
        gnu.write('set size square\n')
        gnu.write('set output \'out/plot4b.eps\'\n')
        gnu.write('set xlabel \'(G/Gmax)^2\'\n')
        gnu.write('set logscale y\n')
        gnu.write('plot ')
        for i in range(len(self.d20)):
            gnu.write('\'out/test.out\' i %i u 1:2:3 ti \'%s\' w err lt %i,\'out/sim.out\' i %i u 1:3 noti w li  lt %i' %(i,self.d20[i],i+1,i,i+1))
            if(i!=len(self.d20)-1):
                gnu.write(',')
        gnu.write('\n')

        gnu.close()
        os.system('gnuplot gnu/gnu.gp')
        
        os.system('arraygraph.py 3 8 0 0 0 0 `ls '+self.datapath+'/raw*.eps` '+self.datapath+'/diff.eps out/plot0.eps out/plot1.eps out/plot1b.eps out/plot2.eps out/plot3.eps out/plot3b.eps out/plot4.eps out/plot4b.eps')
        os.system('mv summary.pdf analysis.pdf')

