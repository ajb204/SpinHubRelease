#123
#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3

###################################
#Lets solve the assignment problem.
#5th April 2019
import enum
import os,sys,numpy,copy,re
import assign.assFrame
import nmrglue as ng
from decon.Frames.deconFrame import readpeaklist,findnear_index,connEntry
import matplotlib
from SettingsUnidec import Parse

#main peak entry class
#needs to map to connEntry class from decon
class peakEntry():
    def __init__(self,test):
        self.name=test[0]
        self.f1=float(test[1])
        self.f2=float(test[2])
        self.f3=float(test[3])
        self.f3p=float(test[3])
        self.inty=float(test[4])

        if(len(self.name.split('_'))==1):
            self.pk=self.name
        else:
            self.pk=self.name.split('_')[0]
            self.ind=self.name.split('_')[1]
        self.tp=''




#spectrum class: contains raw data, peak information in both 2d and 3d
class spectrum():
    def __init__(self,tp,indir,deconFil,peakFil,ref,f1180):

        self.tp=tp
        self.indir=indir
        self.deconParFile=deconFil
        self.peakFil=peakFil
        try:
            self.ref=float(ref)
        except:
            self.ref=ref
        self.f1180=f1180

        '''if(self.Parse()):       # self.Parse() function is True if input file can be found, otherwise it is False 
            self.MakeInp()      # Make Input Function 
            self.GetPeakList()
            self.GetConn()
            pass'''

        if(self.Parse()):       # self.Parse() function is True if input file can be found, otherwise it is False 
            self.MakeInp()      # Make Input Function 
            self.GetPeakList()  # read 3D peak list. 
            self.GetConn()

            #self.Cluster()
            

    ##################################################################
    #Clustering 3Ds to get 2Ds.

    #------------------------------------------------------#
    #probably need to move this into the unidecproject.
    #get peak positions and names. work out averages
    #find if anyone is nearer to someone else. if so, remap
    #repeat until no-one needs remapping.
    #fac is for 15N divisor for CSP.
    def Cluster(self,fac=10):

        print('Clustering 3D peaks')

        posx,posy,posn=self.GetPositions() #read out xyz and name
        pks=numpy.unique(posn) #get the unique peak names.
        go=0
        while(1==1):
            aveX,aveY,aveN=self.GetAve(pks,posx,posy,posn,fac)  #get average positions based on posn
            close,finish,locca=self.GetMaps(posx,posy,posn,aveX,aveY,aveN,fac) #get remap options
            if(len(close)==0): #if there are no peaks closer to someone else than the mean.
                break
            arga=numpy.argmin(close)  #get the peak who is closest to someone else...
            #print('Moving:',self.peak[locca[arga]].name,'to',finish[arga],'dist:',close[arga])
            posn[locca[arga]]=finish[arga] #remap. go again.
            go+=1
            
        #resort peak list based on any numbers in peak.
        import re
        num=[] #sort by any number in the peak name.
        for pk in aveN:
            num.append(int(re.findall(r'\d+', pk)[0]))
        num=numpy.array(num)
        argo=numpy.argsort(num)
        aveX=aveX[argo]
        aveY=aveY[argo]
        aveN=aveN[argo]

        #write out peaklist for  checking.
        #Also: check peaks not in the 3D list, but in the 2D list.
        #make sure these are added.

        #A RISK: we could move the peaks and confuse the issue, having peaks from 
        #the wrong specturm map to the wrong other spectrum.


        outfile=os.path.join(self.pars['indir'],self.pars['peakfile']+'.max')
        outy=open(outfile,'w')
        for i in range(len(aveX)):
            x=aveX[i];y=aveY[i];n=aveN[i]
            outy.write('%s\t%f\t%f\n' % (n,y,x))  #write name, N,H   
        outy.close()
        print('Adjusted %i peaks. New list created: %s' % (go,outfile))
        return 
    
    #get positions, and do a walk in 3D to find maximum intensity.
    def GetPositions(self):
        posx=[]
        posy=[]
        posn=[]
        #step=[]
        
        #infile=os.path.join(self.pars['indir'],self.peakFil,'.max')
        for ii,pk in enumerate(self.peak):

            self.alias(pk,pk.f3,0)
            self.alias(pk,pk.f2,1)
            self.alias(pk,pk.f1,2)
            k=pk.indexK   #findnear_index(pk.f1,self.index2) #find index2
            j=pk.indexJ   #findnear_index(pk.f2,self.index1) #find index1 
            i=pk.indexI   #findnear_index(pk.f3,self.index0) #find index0

            #inty=self.data[i,j,k]
            #print(pk.name,self.index0[i],self.index1[j],self.index2[k],inty,pk.f3,pk.f2,pk.f1,pk.inty)
            arr,steps=self.CheckMax(i,j,k) #get the local maximum.
            pk.indi=arr[0]
            pk.indj=arr[1]
            pk.indk=arr[2]
            pk.steps=steps
            #if(pk.pk not in indy.keys()):
            #    indy[pk.pk]=[]
            #indy[pk.pk].append(ii)

            loc=self.index2[arr[2]],self.index1[arr[1]]
            posx.append(self.index2[arr[2]])
            posy.append(self.index1[arr[1]])
            posn.append(pk.pk)
            #step.append(steps)
            
        #set maxima.
        posx=numpy.array(posx) 
        posy=numpy.array(posy)
        posn=numpy.array(posn) #peak id with peak entry.
        #step=numpy.array(step)
        return posx,posy,posn

    #get the average x/y for each peak.
    #fac is used for CSP measure.
    def GetAve(self,pks,posx,posy,posn,fac):
        aveX=[]
        aveY=[]
        aveN=[]
        dist=[]
        loc=[]

        for pk in pks: #for each unique peak...
            arr=numpy.where(pk==posn)
            loc.append(arr)
            xx=posx[arr]  #1H
            yy=posy[arr]  #15N
            
            aX=numpy.average(numpy.average(xx))
            aY=numpy.average(numpy.average(yy))
            
            aveX.append(aX)
            aveY.append(aY)
            aveN.append(pk)
            dist.append(((xx-aX)**2.+((yy-aY)/fac)**2.)**0.5)
            #print(pk,aX,aY,step[arr],dist[-1])

        aveX=numpy.array(aveX)  #aveX
        aveY=numpy.array(aveY)  #aveY
        aveN=numpy.array(aveN)  #unique name.
        return aveX,aveY,aveN
       
    #go through all peaks, find ones closer to someone else than their own locus    
    def GetMaps(self,posx,posy,posn,aveX,aveY,aveN,fac):
        #loop again...
        close=[]
        finish=[]
        locca=[]
        for ii in range(len(posx)): #for each peak. Who are we nearest to?
            dx=posx[ii]-aveX
            dy=posy[ii]-aveY

            dr=(dx**2+(dy/fac)**2)
            argy=numpy.argmin(dr) #find nearest...
            if(posn[ii]!=aveN[argy]): #if this maps to someone different...
                #argo=numpy.where(aveN==posn[ii])  #flag t
                #dr[argo]**0.5 #self-distance, larger than new distance.
                #print(self.peak[ii].name,posn[ii],aveN[argy],dr[argy]**0.5)
                close.append(dr[argy])    #closest distnace...
                finish.append(aveN[argy]) #map to this peak.
                locca.append(ii)          #location of peak that is going to get mapped.
        return numpy.array(close),finish,locca


    #for given i and dimension, get the triplet i-1,i,i+1
    #if less than zero or above max, adjust.
    def GetTriplet(self,i,dim):
        iv=i-1,i,i+1
        iv=numpy.array(iv,dtype=numpy.int64)
        m1=iv<0
        iv[m1]+=self.specsize[dim] #roll by data size.
        m2=iv>=self.specsize[dim]
        iv[m2]-=self.specsize[dim]
        return iv

    #we will need to use these.
    ind3x3={}
    for ii in range(27):
        ind3x3[ii]=ii//9,ii%9//3,ii%3
    one3=1,1,1
    one3=numpy.array(one3)

    #for the given place, get the 9 steps around given ijk
    #find the max. If the max is in the centre, we are done!
    def Step3D(self,arr):
        #print('y',arr)
        iv=self.GetTriplet(arr[0],0)
        jv=self.GetTriplet(arr[1],1)
        kv=self.GetTriplet(arr[2],2)
        #inty=self.data[iv,jv,kv] #doesn't work.
        inty=self.data[iv,:,:][:,jv,:][:,:,kv]  #not pretty? surely can do this better?
        
        #print(numpy.fabs(inty).shape)
        argy=numpy.argmax(numpy.fabs(inty)) #get max, but argmax unravels.
        ind=self.ind3x3[argy] #map the unraveled max to a 3x3 index. 
        
        if((ind==self.one3).all()): #if central spot is max...
            return arr,True
        arr=numpy.array((iv[ind[0]],jv[ind[1]],kv[ind[2]])) #update.
        return arr,False

        #verify the two are the same.
        #print(numpy.max(numpy.fabs(inty)))
        #print(self.data[iv[ix],jv[jx],kv[kx]])
        #return ind
    
    def CheckMax(self,i,j,k):
        arr=(i,j,k)
        arr=numpy.array((i,j,k))
        #print('Starting intensity:',self.data[arr[0],arr[1],arr[2]])
        steps=0
        while(1==1):
            steps+=1
            #print(arr)
            arr,maxy=self.Step3D(arr)
            if(maxy):
                #print ("Success!")
                break
        #
        #print('Final intensity:',self.data[arr[0],arr[1],arr[2]])
        return arr,steps
    
    #END reclustering 3Ds to find 2D
    ##################################################################



    def GetPeakList(self):  #read a 3D peak list, save peak entries in peak.
        #if('peakfile' not in self.pars.keys()):
        #    print('Cannot find peakfile')
        #    return
        if('indir' not in self.pars.keys()):
            print('Cannot find indir')
            return
        #infile=os.path.join(self.pars['indir'],self.pars['peakfile'])
        infile=os.path.join(self.pars['indir'],self.peakFil)


        if(os.path.exists(infile)==0):
            print('File does not exist:',infile)
            return
        inny=open(infile)
        self.peak=[]
        
        for line in inny.readlines():
            test=line.split()
            self.peak.append(peakEntry(test))   # peak entry class defines the f1, f2, f3 and f3p (aliased 3rd frequency) for each test variable

         
                
        
        print('Read in:',len(self.peak),'3D cross peaks')



    def GetConn(self):
        infile=os.path.join(self.pars['indir'],self.peakFil)
        inny=open(infile)
        self.dim=3
        self.conn_data=[]
        for line in inny.readlines():
            test=line.split()
            self.conn_data.append(connEntry(test,sym='n',peak=self.peak2D,dim=self.dim))
        inny.close()


        self.noeTags=[]
        for i in range(len(self.conn_data)):
            self.noeTags.append(self.conn_data[i].tag)


    def CheckPath(self):
        #print(self.pars['indir'])
        #print('a',self.pars['indir'].split("/"))
        tast=self.pars['indir'].split("/")
        for i in range(len(tast)):
            ii=len(tast)-i-1
            test=os.path.join(os.getcwd(),self.pars['indir'].split("/")[ii],self.pars['infile'])
            if(os.path.exists(test)==1):
                self.pars['indir']=os.path.join(os.getcwd(),self.pars['indir'].split("/")[ii])
                return test
        print('Cannot find directory',self.pars['indir'])
        sys.exit(100)

    def GetProjFile(self,projdir):
        projname=self.labb[-1]+'.'+self.labb[-2]+'.dat'
        projfile=os.path.join(projdir,'projections',projname)

        print(self.labb)
        print(projname,projfile)
        
        if(os.path.exists(projfile)):
            return projfile
        
        projname=self.labb[-2]+'.'+self.labb[-1]+'.dat'
        projfile=os.path.join(projdir,'projections',projname)
        if(os.path.exists(projfile)):
            return projfile

        print(self.labb)
        print(projname,projfile)
        #if(self.tp=='hncocanh' or self.tp=='hncanh'):
        #    projfile=os.path.join(projdir,'projections/'+'Hx'+'.'+'Nx'+'.dat')
        #else:
        #    projfile=os.path.join(projdir,'projections/'+'H'+'.'+'N'+'.dat')
        #try:
        #    self.dic_proj,self.data_proj=ng.pipe.read(projfile)
        #except:
        #    projfile=os.path.join(projdir,'projections/'+'H1'+'.'+'N15'+'.dat')
        #    self.dic_proj,self.data_proj=ng.pipe.read(projfile)

        
    def MakeInp(self):
        main_dir = self.pars['indir']
        if main_dir == './':
             main_dir = os.getcwd()+self.indir[1:]
             self.pars['indir'] = main_dir
        infile=os.path.join(main_dir,self.pars['infile'])
        if(os.path.exists(main_dir)==0):
            infile=self.CheckPath()
            #print(os.path.exists(infile))
        #infile=os.path.join(self.pars['indir'],'out/MyList.out')
        self.dic,self.data=ng.pipe.read(infile)         # self.data is a 3D matrix of intensities, self.dic contains variables
        uc0 = ng.pipe.make_uc(self.dic,self.data,dim=0)
        uc1 = ng.pipe.make_uc(self.dic,self.data,dim=1)
        uc2 = ng.pipe.make_uc(self.dic,self.data,dim=2)
        ord=self.dic['FDDIMORDER']
        lab1=self.dic['FDF1LABEL']
        lab2=self.dic['FDF2LABEL']
        lab3=self.dic['FDF3LABEL']
        lab=lab1,lab2,lab3
        self.labb=lab[int(ord[2])-1],lab[int(ord[1])-1],lab[int(ord[0])-1]
        self.dmax=numpy.max(self.data)  #get max intensity in spectrum

        Size=self.data.shape
        self.noise=self.dmax*float(self.pars['thresh']  )
        self.specsize=Size
        self.uc0max=uc0.ppm(0)
        self.uc0min=uc0.ppm(Size[0]-1)
        self.uc1max=uc1.ppm(0)
        self.uc1min=uc1.ppm(Size[1]-1)
        self.uc2max=uc2.ppm(0)
        self.uc2min=uc2.ppm(Size[2]-1)
        print('--------------------------------------------------')
        print('Reading:',infile)
        print("Spectrum dimensions (pts): ",Size)   #print(the spectral dimensions)
        print("Labels: ",self.labb)
        print("dimension 0 limits (ppm): ", self.uc0min, self.uc0max)  #carbon)
        print("dimension 1 limits (ppm): ", self.uc1min, self.uc1max)  #direct)
        print("dimension 2 limits (ppm): ", self.uc2min, self.uc2max)  #direct)
        #print('Maximum Intensity:',self.dmax)
        self.index0=[]#make index of carbon chemical shifts for index 0
        for i in range((Size[0])):
            self.index0.append((uc0.ppm(0)-i*(-uc0.ppm(Size[0]-1)+uc0.ppm(0))/(Size[0]-1)))
        self.index1=[]#make index of carbon chemical shifts for index 1
        for i in range((Size[1])):
            self.index1.append((uc1.ppm(0)-i*(-uc1.ppm(Size[1]-1)+uc1.ppm(0))/(Size[1]-1)))
        self.index2=[]#make index of carbon chemical shifts for index 2
        for i in range((Size[2])):
            self.index2.append((uc2.ppm(0)-i*(-uc2.ppm(Size[2]-1)+uc2.ppm(0))/(Size[2]-1)))
        self.index0=numpy.array(self.index0)
        self.index1=numpy.array(self.index1)
        self.index2=numpy.array(self.index2)
        self.alias0=numpy.max((self.uc0min,self.uc0max))-numpy.min((self.uc0min,self.uc0max))+numpy.fabs(self.index0[0]-self.index0[1])
        self.alias1=numpy.max((self.uc1min,self.uc1max))-numpy.min((self.uc1min,self.uc1max))+numpy.fabs(self.index1[0]-self.index1[1])
        self.alias2=numpy.max((self.uc2min,self.uc2max))-numpy.min((self.uc2min,self.uc2max))+numpy.fabs(self.index2[0]-self.index2[1])
        self.YY,self.XX,self.ZZ=numpy.meshgrid(self.index1,self.index0,self.index2)
        ########GET PROJECTION##############
        tmp=infile.split('/')[:-1]
        projdir='/'
        for tm in tmp:
            if(tm!='.'):
                projdir=os.path.join(projdir,tm)
        ##projfile=os.path.join(projdir,'projections/'+lab2+'.'+lab3+'.dat')
        ##TEMPORARY FIX OF BUG - AT THE MOMENT IT IS USING THE CH PLANE NOT NH IN SPECTRA THAT NEED CAREFUL TRANSPOSING

        #print(self.labb)
        projfile=self.GetProjFile(projdir)
        self.dic_proj,self.data_proj=ng.pipe.read(projfile)

        uc0_proj = ng.pipe.make_uc(self.dic_proj,self.data_proj,dim=0)
        uc1_proj = ng.pipe.make_uc(self.dic_proj,self.data_proj,dim=1)
        ord_proj=self.dic['FDDIMORDER']
        lab1_proj=self.dic['FDF1LABEL']
        lab2_proj=self.dic['FDF2LABEL']
        lab3_proj=self.dic['FDF3LABEL']
        lab_proj=lab1_proj,lab2_proj,lab3_proj
        #print(ord_proj)
        #print(lab_proj)
        self.labb_proj=lab_proj[int(ord_proj[1])-1],lab_proj[int(ord_proj[0])-1]
        #self.dmax=numpy.max(self.data)  #get max intensity in spectrum
        #print('dmax:',self.dmax)
        Size_proj=self.data_proj.shape
        #self.noise=self.dmax*float(self.pars['thresh']  )
        self.specsize_proj=Size
        self.uc0max_proj=uc0_proj.ppm(0)
        self.uc0min_proj=uc0_proj.ppm(Size_proj[0]-1)
        self.uc1max_proj=uc1_proj.ppm(0)
        self.uc1min_proj=uc1_proj.ppm(Size_proj[1]-1)
        print("Projection dimensions (pts): ",Size_proj)   #print(the spectral dimensions)
        print("Labels: ",self.labb_proj)
        print("dimension 0 limits (ppm): ", self.uc0min_proj, self.uc0max_proj)  #carbon)
        print("dimension 1 limits (ppm): ", self.uc1min_proj, self.uc1max_proj)  #direct)
        #print('Maximum Intensity:',self.dmax)
        self.index0_proj=[]#make index of carbon chemical shifts for index 0
        for i in range((Size_proj[0])):
            self.index0_proj.append((uc0_proj.ppm(0)-i*(-uc0_proj.ppm(Size_proj[0]-1)+uc0_proj.ppm(0))/(Size_proj[0]-1)))
        self.index1_proj=[]#make index of carbon chemical shifts for index 1
        for i in range((Size_proj[1])):
            self.index1_proj.append((uc1_proj.ppm(0)-i*(-uc1_proj.ppm(Size_proj[1]-1)+uc1_proj.ppm(0))/(Size_proj[1]-1)))
        self.index0_proj=numpy.array(self.index0_proj)
        self.index1_proj=numpy.array(self.index1_proj)
        self.YY_proj,self.XX_proj=numpy.meshgrid(self.index1_proj,self.index0_proj)
        #print(self.XX_proj.shape,self.YY_proj.shape,self.data_proj.shape)
        ######NOW DO PEAKS#####
        peakListLocation=os.path.join(self.pars['indir'],self.pars['peakfile'])
        if(peakListLocation.split('.')[-1]!= 'list'): #reset 2D peak file name if not of type list.
            peakListLocation = peakListLocation.split('.')[0]+'raw/test.ft3.list'

        """
        outp = open('error.txt', 'w')
        outp.write('%s \n' % self.pars['indir'])
        outp.write('%s \n' % self.pars['peakfile'])
        outp.write('%s' % peakListLocation)
        outp.close()
        """

        self.peak2D=readpeaklist(peakListLocation)
        for p in range(len(self.peak2D)): #take 2D peak list, and alias all 3 to be in range.
            self.alias(self.peak2D[p],self.peak2D[p].y,0)
            self.alias(self.peak2D[p],self.peak2D[p].y,1)
            self.alias(self.peak2D[p],self.peak2D[p].x,2)
        self.PEAK=1
        self.pkIdx=[] #index of peak positions
        self.pkSlice1D=[] #1D slices
        for p in range(len(self.peak2D)):
            ptC=self.peak2D[p].indexJ #carbon
            ptH=self.peak2D[p].indexK #proton
            self.pkIdx.append((ptC,ptH))
            self.pkSlice1D.append(self.data[:,ptC,ptH])

    #for a given dimension, find relevant alias information, 
    #max, min, sw, such that fold will be (max-min)+dd
    def GetAliasRange(self,dim):
        if(dim==0):
            dd=numpy.fabs(self.index0[1]-self.index0[0] )
            dmax=self.uc0max
            dmin=self.uc0min
            vals=self.index0
            return dd,dmax,dmin,vals
        if(dim==1):
            dd=numpy.fabs(self.index1[1]-self.index1[0] )
            dmax=self.uc1max
            dmin=self.uc1min
            vals=self.index1
            return dd,dmax,dmin,vals
        if(dim==2):
            dd=numpy.fabs(self.index2[1]-self.index2[0] )
            dmax=self.uc2max
            dmin=self.uc2min
            vals=self.index2
            return dd,dmax,dmin,vals
        if(dim==3):
            dd=numpy.fabs(self.index3[1]-self.index3[0])
            dmax=self.uc3max
            dmin=self.uc3min
            vals=self.index3
            return dd,dmax,dmin,vals
        #whoops!
        return False,False,False,False
        


    #for given ppm for specified peak in dimension dim
    #fold until we are in the range of the spectrum.
    def alias(self,peak,ppm,dim):
        
        dd,dmax,dmin,vals=self.GetAliasRange(dim) #for specified dimension, get alias values
        if(dd==False):
            return -1

        #print(ppm)
        while(ppm>=dmax):
            ppm-=(dmax-dmin+dd)
        while(ppm<=dmin):
            ppm+=(dmax-dmin+dd)
        
        #print(ppm)
        i=findnear_index(ppm,vals)

        if(dim==0):
            peak.indexI=i
            peak.ppmI=vals[i]
        elif(dim==1):
            peak.indexJ=i
            peak.ppmJ=vals[i]
        elif(dim==2):
            peak.indexK=i
            peak.ppmK=vals[i]
        elif(dim==3):
            peak.indexL=i
            peak.ppmL=vals[i]

    def Parse(self):
        self.pars={}
        infile=os.path.join(self.indir,self.deconParFile)
        if self.indir == './':
            self.indir = os.getcwd()
        if(os.path.exists(infile)==0):
            print('Cannot find:',infile)
            return 0
        inny=open(os.path.join(self.indir,self.deconParFile))
        for line in inny.readlines():
            test=line.split()
            if(len(test)==3):
                if(test[1]=='='):
                    self.pars[test[0]]=test[2]

        return 1
    
    def findnear_index(self,test,array):
        #array = numpy.asarray(array)
        idx = (numpy.abs(array - test)).argmin()
        return idx
        #return array[idx]

    #spectrum, peak name, specify target ppm and return intensity
    def GetIntensity(self,pk,targ):      
            pkl=-1
            for i,n in enumerate(self.peak2D):
                    #print(pk,n.name)
                    if(n.name==pk):
                        pkl=i
                        break
            if(pkl==-1):
                    print('peak not found')
                    return 0,0
            
            #print(pkl)
            ptC=self.pkIdx[pkl][0]
            ptH=self.pkIdx[pkl][1]
            Xs=self.XX[:,ptC,ptH]
            Zs=self.data[:,ptC,ptH] #extract the relevant 2d slice
            
            ptY=self.findnear_index(targ,Xs-self.ref)
            #print(Xs,Zs)
            #print(Zs.shape)
            #print(ptY)
            #print(TARG,Xs[ptY],Zs[ptY])
            #print('pty',ptY)
            return Xs[ptY],Zs[ptY],Zs[ptY]/self.noise
    
    def GetSlice(self,pk):      
            pkl=-1
            for i,n in enumerate(self.peak2D):
                    #print(pk,n.name)
                    if(n.name==pk):
                        pkl=i
                        break
            if(pkl==-1):
                    print('peak not found')
                    return 0,0
            
            #print(pkl)
            ptC=self.pkIdx[pkl][0]
            ptH=self.pkIdx[pkl][1]
            Xs=self.XX[:,ptC,ptH]
            Zs=self.data[:,ptC,ptH] #extract the relevant 2d slice
            #TARG=58
            #ptY=findnear_index(targ,Xs-mol.spec[lab].ref)
            #print(Xs,Zs)
            #print(Zs.shape)
            #print(ptY)
            #print(TARG,Xs[ptY],Zs[ptY])
            #print('pty',ptY)
            return Xs-self.ref,Zs      
    
    def FindPeakMaxima(self,peak,width,h):
        X,Z=self.GetSlice(peak)   #extract slice from raw data
        dx=numpy.fabs(X[1]-X[0])  #work out spacing.
        dist=width/dx #1ppm divided by dwell space is minimum distance
        import scipy
        pks,prop=scipy.signal.find_peaks(Z,distance=dist,width=(None,None),height=self.noise*h)
        return X[pks],Z[pks]
    

    #should be able to do this without specifying peaks.
    def GetSets(self,peaks,cspH=0.1,cspN=1.,verb=True):
        #SPECTRUMCLASS
        
        
        loc=[]
        pks=[]
        for pk in peaks: #go over user inputted list...
            
            #self.peak.append(peakEntry(test))
            tig=0
            for pok in self.peak:
                pik=pok.name.split('_')[0]
                if(pik==pk):
                    tig=1
                    break
            if(tig==0):
                continue
            
            x=pok.f1 #1H
            y=pok.f2 #15N
            #print(x,y)
            #pkl=-1   #match this to someone in peak2D
            #for i,n in enumerate(self.peak2D):
            #        #print(pk,n.name)
            #        if(n.name==pk):
            #            pkl=i
            #            break
            #if(pkl==-1):
            #    continue

            #x=self.peak2D[pkl].x #save x
            #y=self.peak2D[pkl].y #save y
            l=(x,y)
            l=numpy.array(l)
            loc.append(l)
            pks.append(pk)
        loc=numpy.array(loc)
        pks=numpy.array(pks)

        sets=[] #now cluster.
        for i,pk in enumerate(pks): #from each position to all other positions
            dH=((loc[i,0]-loc[:,0]))
            dN=((loc[i,1]-loc[:,1]))
            r2=dH**2+dN**2.  #square distanace from point to all other points

            ang=numpy.arctan2(dN,dH)
            rLim2=(cspH*numpy.cos(ang))**2.+ (cspN*numpy.sin(ang))**2.
            mask=(r2<rLim2)*(numpy.fabs(dH)<cspH)*(numpy.fabs(dN)<cspN)

            #condition for ellipse intersection
          
            if(numpy.sum(mask*1.)>1): #if there is more than 1 (1 being itself)...
                toadd=pks[mask] #get a list of the affected peaks...
                #r2s=r2[mask]
                #r2l=rLim2[mask]
                #hh=dH[mask]
                #nn=dN[mask]
                #for (t,h,n,r,l) in zip(toadd,hh,nn,r2s,r2l):
                #    print(pk,t,h,n,r**0.5,l**0.5)

                jj=-1
                for a in toadd:  #for each member...
                    for j in range(len(sets)): #am I already in any sets?
                        if(a in sets[j]): #yes I am!
                            jj=copy.deepcopy(j) #set value to that of j.
                            break 
                if(jj==-1): #none of the peaks are in a pre-existing set.
                    line=sorted(toadd)  #make a new set.
                    #print('adding:',toadd)
                    sets.append(line)
                else: #at least one of the peaks is in a set, add all to the set, if not there.
                    for a in toadd:
                        if(a not in sets[jj]):
                            #print('appending:',a)
                            sets[jj].append(a)
        
        if(verb):
            for i,s in enumerate(sets): #tell everyone what wonderful sets we have.
                print ('Set ',i+1,':',s)
            
        return sets




class molecule():
    def __init__(self):
        self.LOADASSIGN=0  #have assignments been loaded?

        self.tolHNCO=0.02*3   #tolerance for hnco
        self.tolHNCA=0.2*2   #tolernace for hnca
        self.tolHNCO=0.3   #tolerance for hnco
        self.tolHNCA=0.3   #tolernace for hnca
        self.tolHNCACB=0.3 #tolerance for hncacb
        self.tolHNCANH=0.3  # tolerence for hncanh

        #self.Hfudge=(4.77-4.27)
        self.Hfudge= 0
        #chemical shift discriminators for residue classificaiton
        #against BMRB
        self.tolSin=5.  #factor by which we have awinner
        self.tolSin=20.  #factor by which we have awinner
        self.tolMax=0.6 #0.6
        self.tolMin=0.05 #0.1

        self.tolMatch=0.2*2. #tolerance for matching edges - exclude above this value


        self.tolMatch=0.3




        #self.refSpec='hncaco'

        #self.template='4mjh.pdb.cs'
        #self.template='hsp27.bmrb'
        #self.template='1UBQ.pdb.cs'

        self.refSpec='hnco'
        self.template=''



        self.peak={}     #will contain peak library
        self.G1edges={}  #will contain graphs
        self.Optedges={} #will contain possible peak matches
        self.spec={}     #will contain spectra

        self.HNCOmed=0   #median values (hnco-hncaco)
        self.HNCAmed=0   #median values (hnca-hncoca)
        self.HNCAHNCACBmed=0   #median values (hnca-hncoca)
        self.HNCACBmed=0 #median values (hncacb-hncocacb)
        self.HNCANHmed=0
        self.CTOCSYmed=0
        self.HTOCSYmed=0

        self.p1to3 = {
            'A': 'ALA', 'C': 'CYS', 'D': 'ASP',
            'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS',
            'I': 'ILE', 'K': 'LYS', 'L': 'LEU', 'M': 'MET',
            'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
            'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP',
            'Y': 'TYR',
        } #translate 1name to 3name

        self.p3to1 = dict((x[1], x[0]) for x in self.p1to3.items()) #translate 3name to 1name

        self.SetBMRB()

    def GetSeq(self,infile):
        self.seq={}
        inny=open(infile)

        if(infile=='sequence.txt'):
            cnt=1
            for line in inny.readlines():
                test=list(line.split('\n')[0])
                self.FirstResidue=1
                for i in range(len(test)):
                    resi = cnt
                    resn = test[i]
                    self.seq[resi]=resn
                    cnt+=1
                    

        #cnt=0
        elif(infile[-4:]=='.seq'):
            cnt=1
            self.FirstResidue=1  #default to first residue is 1
            for line in inny.readlines():

                test=line.split()
                if(len(test)>0):
                    if(test[0]=='start:' or test[0]=='Start:'):
                        cnt=int(test[1])  #adjust first residue 
                        self.FirstResidue=int(test[1])
                    else:
                        for te in test:
                            for t in te:
                                # print(t)
                                resi=cnt
                                resn=t
                                #print(resi,resn,self.p3to1[resn])
                                self.seq[resi]=resn
                                cnt+=1
                                #cnt+=1
                    #del self.seq[1] #remove first residue


        else:
            for line in inny.readlines():
                test=line.split()
                if(len(test)>0):
                    #if(cnt>0):
                    resi=int(test[4])
                    resn=test[6]
                    #print(resi,resn,self.p3to1[resn])
                    self.seq[resi]=self.p3to1[resn]
                    #cnt+=1
                    #del self.seq[1] #remove first residue
            self.FirstResidue=numpy.min(list(self.seq.keys()))

        
        self.resi=sorted(self.seq.keys()) #residue numbers, sorted.


        self.G1_nodes=[]
        for i in range(len(self.resi)): #residue names
            self.G1_nodes.append((str(self.resi[i])+self.seq[self.resi[i]]))
            #print(self.G1_nodes[-1])
        #sys.exit(100)
        self.G1_noes={}
        for i in range(len(self.resi)):
            #self.G1_nodes.append((self.resi+self.seq[self.resi[i]]))
            if(self.seq[self.resi[i]]!='P'):
                self.G1_noes[self.G1_nodes[i]]=[]
                if(i!=len(self.resi)-1 and self.seq[self.resi[i+1]]!='P'):
                    self.G1_noes[self.G1_nodes[i]].append( (self.G1_nodes[i+1],1,'f'))
                if(i!=0  and self.seq[self.resi[i-1]]!='P'):
                    self.G1_noes[self.G1_nodes[i]].append( (self.G1_nodes[i-1],1,'b'))
        print('Number of residues:',len(self.resi))
        print('First residue:',self.FirstResidue)

    #add a peak list from a spectrum to the peak library
    def AddSpec(self,spec):
        self.spec[spec.tp]=spec #add new spectrum instance.
        for pk in spec.peak: #now copy in the 3D peaks from the spectrum into molecule list.
            if(pk.pk not in self.peak.keys()): #if peakname not already there...
                self.peak[pk.pk]={}
            if(spec.tp not in self.peak[pk.pk].keys()):  #if spectrum type not there...
                self.peak[pk.pk][spec.tp]=[]
            self.peak[pk.pk][spec.tp].append(pk)  #append peak entry.



    #print(to screen the peak library)
    def listy(self):
        keys=self.peak.keys()
        keys=sorted(keys)
        for pk in keys:
            for spec in self.peak[pk]:
                for i,pk3 in enumerate(self.peak[pk][spec]):
                    print(pk,spec,pk3.name,pk3.f1,pk3.f2,pk3.f3,pk3.tp,pk3.inty)

    #go through peak and rename it all, save old name
    def RenamePeaks(self):
        for pk in self.peak.keys():
            for spec in self.peak[pk]:
                for pk3 in self.peak[pk][spec]:
                    lab=self.renameAss(spec,pk,pk3.name,pk3.tp)
                    pk3.oldname=pk3.name
                    pk3.name=lab

    #rename based on assignment
    def renameAss(self,spec,pk,v1,tp):
        if(pk not in self.results.keys()):
            return v1
        if(len(self.results[pk])!=1):
            return v1
        # print(pk)
        # print(self.results[pk])
        ass=self.results[pk][0]
        resi=int(re.findall(r'\d+',ass)[0])
        resn=ass.split(str(resi))[1]
        resn3=self.p1to3[resn]
        lead='%s%iHN-N-'% (resn3,resi)
        if(spec=='hnco'):
           return '%s%i%s' % (lead,resi-1,'CO')
        elif(spec=='hncaco'):
            if(tp=='main'):
                return '%s%i%s' % (lead,resi,'CO')
            else:
                return '%s%i%s' % (lead,resi-1,'CO')
        elif(spec=='hnca'):
            if(tp=='main'):
                return '%s%i%s' % (lead,resi,'CA')
            else:
                return '%s%i%s' % (lead,resi-1,'CA')
        elif(spec=='hncoca'):
            #if(tp=='main'):
            return '%s%i%s' % (lead,resi-1,'CA')
        elif(spec=='hncacb'):
            if(tp=='PosMax'):
                return '%s%i%s' % (lead,resi,'CA')
            elif(tp=='NegMax'):
                return '%s%i%s' % (lead,resi,'CB')
            elif(tp=='PosMin'):
                return '%s%i%s' % (lead,resi-1,'CA')
            elif(tp=='NegMin'):
                return '%s%i%s' % (lead,resi-1,'CB')
        elif(spec=='hncocacb'):
            if(tp=='NegMin'):
                return '%s%i%s' % (lead,resi-1,'CB')
        return v1

    def SpecAdjust(self,spec):
        return self.spec[spec].ref
        """
        if(spec=='hncaco'):
            return self.HNCOmed
        if(spec=='hncacb'):
            return self.HNCAHNCACBmed
        if(spec=='hncoca'):
            return self.HNCAmed
        if(spec=='hncocacb'):
            return self.HNCACBmed
        if(spec=='cbcaconh' or spec=='hncocacb'):           ## CB
            return self.HNCACBmed          ## CB
        if(spec=='ctocsy'):
            return self.CTOCSYmed
        if(spec=='hcconh'):
            return self.HTOCSYmed
        
        return 0
        """

    def save_peaks(self,outfile):
        keys=self.peak.keys()
        keys=sorted(keys)
        outy=open(outfile,'w')
        cnt=0
        for pk in keys:
            for spec in self.peak[pk]:
                for i,pk3 in enumerate(self.peak[pk][spec]):
                    if(self.LOADASSIGN):
                        outy.write('%s %s %s %f %f %f %e %s %s\n' % (spec,pk,pk3.name,pk3.f1,pk3.f2,pk3.f3,pk3.inty,pk3.tp,pk3.oldname))
                    else:
                        outy.write('%s %s %s %f %f %f %e %s\n' % (spec,pk,pk3.name,pk3.f1,pk3.f2,pk3.f3,pk3.inty,pk3.tp))
                    cnt+=1
        outy.close()
        print('Saved',cnt,'peaks in',outfile)

        for spec in self.spec.keys():
            for pk in self.peak.keys():
                if(spec in self.peak[pk].keys()):
                    for pk3 in self.peak[pk][spec]:
                        test=pk3.name
                        tig=0;
                        for cn in self.spec[spec].conn_data:
                            if(cn.p1==test):
                                tig=1
                        if(tig==0):
                            print('PEAK NOT IN ORIGINAL FILE')
                            print(pk3.name)
                            print(pk3.name,pk3.f1,pk3.f2,pk3.f3,pk3.inty)
                            f3=pk3.f3 #adjust peak position by offset

                            f3+=self.SpecAdjust(spec) #adjust according to medians
       
                            test=(pk3.name,pk3.f1,pk3.f2,f3,pk3.inty) #add peak
                            self.spec[spec].conn_data.append(connEntry(test,sym='n',peak=self.spec[spec].peak2D,dim=self.spec[spec].dim))


            remove=[]
            for i,cn in enumerate(self.spec[spec].conn_data):
                tig=0
                for pk in self.peak.keys():
                    if(spec in self.peak[pk].keys()):
                        for pk3 in self.peak[pk][spec]:
                            if(cn.p1==pk3.name): #update connectivity with peak list
                                tig=1
                                cn.f1=pk3.f1
                                cn.f2=pk3.f2
                                cn.tp=pk3.tp
                                cn.s1=pk3.inty

                                cn.f3=pk3.f3+self.SpecAdjust(spec) #adjust according to medians

                                #if(spec=='hncaco'): #correct transfer by the offset
                                #    cn.f3=pk3.f3+self.HNCOmed
                                #if(spec=='hncoca'):
                                #    cn.f3=pk3.f3+self.HNCAmed
                                #if(spec=='hncocacb'):
                                #    cn.f3=pk3.f3+self.HNCACBmed
                                #if(spec=='cbcaconh'):                   ## CB
                                #    cn.f3=pk3.f3+self.HNCACBmed         ## CB
                                cn.s1=pk3.inty #unpdate signal intensity
                if(tig==0):
                    remove.append(i)
            remove=sorted(remove,reverse=True)
            for rem in remove:
                print('removing ',self.spec[spec].conn_data[rem].p1)
                self.spec[spec].conn_data.pop(rem)


            infile=os.path.join(self.spec[spec].pars['indir'],self.spec[spec].peakFil)
            print('Copying ',infile,' to make ',infile+'.back')
            os.system('cp '+infile+' '+infile+'.back')

            print('Overwriting',infile)
            outy=open(infile,'w')
            for cn in self.spec[spec].conn_data:
                outy.write('%s\t%f\t%f\t%f\t%e\n' % (cn.p1,cn.f1,cn.f2,cn.f3,cn.s1))
            outy.close()


            if(self.LOADASSIGN): #if assignments are loaded, then rename
                outy=open(infile+'.assign','w')
                for cn in self.spec[spec].conn_data:
                    pk=cn.v1 #peak label
                    lab=self.renameAss(spec,pk,cn.v1,cn.tp)
                    outy.write('%s\t%f\t%f\t%f\t%e\n' % (lab,cn.f1,cn.f2,cn.f3,cn.s1))
                outy.close()

            #print(infile)
            #print(os.path.exists(infile))
            #outy=open(infile,'a')
            #outy.write('%s\t%f\t%f\t%f\t%e\n' % (pk3.name,pk3.f1,pk3.f2,pk3.f3,pk3.inty))
            #outy.close()
            #add it to original list



    #load a peak list
    def load_peaks(self,infile):

        inny=open(infile)
        self.peak={}
        cnt=0
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                if(len(test)>=7):
                    spec=test[0]
                    pk=test[1]
                    inst=peakEntry(test[2:7])
                    if(len(test)==8):
                        inst.tp=test[7]
                        if(spec=='hncacb' or spec=='hncocacb'):
                            if(inst.tp not in ('PosMin','PosMax','NegMin','NegMax')):
                               print('ERROR READING LINE - UNRECOGNISED LABEL:')
                               print('should be PosMin, PosMax, NegMin or NegMax')
                               print(test)
                        else:
                            if(inst.tp not in ('main',)):
                               print('ERROR READING LINE - UNRECOGNISED LABEL:')
                               print('should be either main or nothing')
                               print(test)

                    if(pk not in self.peak.keys()):
                        self.peak[pk]={}
                    if(spec not in self.peak[pk].keys()):
                        self.peak[pk][spec]=[]
                    self.peak[pk][spec].append(inst)
                    cnt+=1

        print('Read in ',cnt,'peaks from',infile)
        self.assSpec()    #link up the peaks
        self.EdgeScreen() #screen to create graphs
        
    def AddManualPeaks(self):
        if('add' not in self.manual.keys()):
            return
        for test in self.manual['add']:
                if(len(test)!=4):
                    print('Wrong number of entries. cannot add.')
                    print(test)
                    print('We need spec,peak,f3,inty')
                    continue
                spec=test[0]
                pk=test[1]
                f3=test[2]
                inty=test[3]
                if(pk in self.peak.keys()):
                    if(spec in self.peak[pk].keys() and len(self.peak[pk][spec])>0):
                        f1=self.peak[pk][spec][0].f1
                        f2=self.peak[pk][spec][0].f2
                        name=pk+'_A'
                        print('Adding:',name,f1,f2,f3,inty,'to',pk,'in',spec)
                        self.peak[pk][spec].append(peakEntry((name,f1,f2,f3,inty)))
                    else:
                        self.peak[pk][spec]=[]
                        go=0
                        for spoc in self.peak[pk].keys():
                            if(len(self.peak[pk][spoc])>0):
                                f1=self.peak[pk][spoc][0].f1
                                f2=self.peak[pk][spoc][0].f2    
                                go=1
                                break
                        if(go==0):
                            print('Trying to add:',name,f1,f2,f3,inty,'to',pk,'in',spec)
                            print('Cannot add peak. Skipping.')
                            continue
                        name=pk+'_A'
                        print('Adding:',name,f1,f2,f3,inty,'to',pk,'in',spec)
                        self.peak[pk][spec].append(peakEntry((name,f1,f2,f3,inty)))

             

    def ReadManual(self):
        self.manual={}
        print('Reading manual adjustments')
        if(os.path.exists('Manual.txt')==False):
            return
        inny=open('Manual.txt')
        for line in inny.readlines():
            test=line.split(',')
            print(test)
            if(len(test)<2):
                continue
            if(test[0][0]=='#'):
                continue
            line=line.replace(' ','')
            print(line)
            test=line.split('\n')[0].split('#')[0].split(',')
            print(test)
            if(test[0] not in self.manual.keys()):
                self.manual[test[0]]=[]
            self.manual[test[0]].append(test[1:])
            print('ManualAdjustment:',test[0],test[1:])
        inny.close()
        print(self.manual)


    #get all peaks from a spectrum, H/N locations
    def GetHNLoc(self,spec):
        loc=[]
        pks=[]
        for pk,specs in self.peak.items():
            if(spec in specs):
                pk3=specs[spec][0]

                x=pk3.f1 #1H
                y=pk3.f2 #15N
                l=(x,y)
                l=numpy.array(l)
                loc.append(l)
                pks.append(pk)
        loc=numpy.array(loc)
        pks=numpy.array(pks)
        return pks,loc

    #should be able to do this without specifying peaks.
    def GetSets(self,spec,verb=True):
        self.SetCSPs()
        pks,loc=self.GetHNLoc(spec)
        #print('    ',spec,self.cspH,self.cspN)
        sets=[] #now cluster.
        for i,pk in enumerate(pks): #from each position to all other positions
            #print('   ',pk)
            dH=((loc[i,0]-loc[:,0]))
            dN=((loc[i,1]-loc[:,1]))
            r2=dH**2+dN**2.  #square distanace from point to all other points

            ang=numpy.arctan2(dN,dH)
            rLim2=(self.cspH*numpy.cos(ang))**2.+ (self.cspN*numpy.sin(ang))**2.
            mask=(r2<rLim2)*(numpy.fabs(dH)<self.cspH)*(numpy.fabs(dN)<self.cspN)

            #condition for ellipse intersection
          
            if(numpy.sum(mask*1.)>1): #if there is more than 1 (1 being itself)...
                toadd=pks[mask] #get a list of the affected peaks...
                #print('toadd',toadd)
                #r2s=r2[mask]
                #r2l=rLim2[mask]
                #hh=dH[mask]
                #nn=dN[mask]
                #for (t,h,n,r,l) in zip(toadd,hh,nn,r2s,r2l):
                #    print(pk,t,h,n,r**0.5,l**0.5)

                tick=0
                jj=-1
                for a in toadd:  #for each member...
                    for j in range(len(sets)): #am I already in any sets?
                        if(a in sets[j]):
                            #print('found',a,'in',sets[j],j)
                            jj=copy.deepcopy(j)
                            tick=1  #yes I am!
                            break 
                #print(jj,tick)
                if(tick==0): #none of the peaks are in a pre-existing set.
                    line=sorted(toadd)  #make a new set.
                    #print('adding:',toadd)
                    sets.append(line)
                    #print('starting:',sets[-1])
                else: #at least one of the peaks is in a set, add all to the set, if not there.
                    for a in toadd:
                        if(a not in sets[jj]):
                            #print (jj)
                            sets[jj].append(a)
                            #print('appending:',a,sets[jj],jj)
        if(verb):
            for i,s in enumerate(sets): #tell everyone what wonderful sets we have.
                print ('Set ',i+1,':',s)
        #sys.exit(10)
        return sets


    #go over the tps. If the name is a number followed by something in aa,
    #save it linked to the appropriate peak.
    def SaveEntries(self,pky,inty,name,tp,ss,aa,spec,stryLine):
                save={} #work out which peakEntries are going to save. back them up
                #(not used)
                for pk in numpy.unique(pky):
                    for a in aa:  #go through each string type to match
                        for i in range(len(inty)):
                            #print(name[i],inty[i],pos[i],tp[i])
                            if(tp[i]==pk+a): #save the is
                                if(pk not in save.keys()):
                                    save[pk]=[]
                                save[pk].append(name[i])
                                break
                    
                pkSave={} #get peak entries together that we're4 going to save, indexed by peak
                for pk,nams in save.items():
                    pkSave[pk]=[]
                    for nam in nams:
                        pok =nam.split('_')[0]
                        for pk3 in self.peak[pok][spec]:
                            if(nam==pk3.name):
                                pkSave[pk].append(copy.deepcopy(pk3))
                                break
                
            
                for pk,pks in pkSave.items():  #go through saved list of peak entries...
                    print('Analysing:',pk)
                    for pk3 in pks: #for each peak...
                        #print('Final:',pk,pk3.name)
                        curr=pk3.name.split('_')[0]
                        if(curr==pk):
                            print('keeping',pk,pk3.name,pk3.f3p,pk3.inty)
                            pass
                            #no action
                        else:
                            #print('kill',pk3.name,'from',curr)
                            print('add',pk3.name,'to',pk,pk3.f3p,pk3.inty)
                            pk4=self.GetPeaksNear(pk,spec,pk3.f3p) #get the matching guy.

                            stry='add,%s,%s,%s,%.2f # Transferred from %s in overlap group %s' % (spec,pk,pk3.f3p,pk4.inty,pk3.name,ss)
                            stryLine.append(stry)   
                        #print(pk,pk3.name)
                    
                    for pk4 in self.peak[pk][spec]: #for each peak currently in there..
                        tick=0
                        for pk3 in pks: #for each peak...
                            if(pk4.name==pk3.name): #
                                tick=1
                        if(tick==1):
                            continue
                        print('scrubbing',pk4.name,'from',pk,pk4.name,pk3.name)
                        stry='remove,%s,%s,%s # Cannot assign in overlap group %s' % (spec,pk,pk4.name,ss)
                        stryLine.append(stry)

    #go over peaks in the overlap group ss
    #populate the lists based on speca.
    #if there is a match of peak position to specb
    #classify that peak as an '(i-1)'
    def AssembleSetInfo(self,ss,speca,specb,thresh,nolabel,blabel,local=False):
            inty=[]
            name=[]
            pos=[]
            pky=[]
            tp=[]
   
            for pk in ss:
                self.ShowPeaks(speca,pk)
                self.ShowPeaks(specb,pk)
                if(speca not in self.peak[pk].keys()):
                    continue
                for pk3 in self.peak[pk][speca]:
                    inty.append(pk3.inty)
                    pos.append(pk3.f3p)
                    name.append(pk3.name)
                    pky.append(pk)

                    tig=0
                    if(local==False):
                        
                        for pok in ss:
                            if(specb not in self.peak[pok]):
                                continue
                            for pk4 in self.peak[pok][specb]:
                                if(numpy.fabs(pk3.f3p-pk4.f3p)<thresh):
                                    tig=1
                                    break
                    else:
                        #for pok in ss:
                        if(specb not in self.peak[pk]):
                            continue
                        for pk4 in self.peak[pk][specb]:
                            if(numpy.fabs(pk3.f3p-pk4.f3p)<thresh):
                                tig=1
                                break


                    if(tig==0):
                        tp.append(nolabel)  #unclassified.
                    else:
                        if(local==False):
                            tp.append(blabel) #classify as a 'bmatch'
                        else:
                            tp.append(pk+blabel) #classify as a 'bmatch'

            for i in range(len(inty)):
                print(name[i],inty[i],pos[i],tp[i])    

            inty=numpy.array(inty)
            pos=numpy.array(pos)
            name=numpy.array(name)
            pky=numpy.array(pky)
            tp=numpy.array(tp)
            argys=numpy.flip(numpy.argsort(numpy.fabs(inty)))
            inty=inty[argys]
            pos=pos[argys]
            name=name[argys]
            pky=pky[argys]
            tp=tp[argys]
            return inty,pos,name,pky,tp

    #for index position ind, look for the label
    #ind+tag (eg 35i). If it's assigned already, skip.
    #if not, assign all positions close in chemical shift
    #(within thresh) to this.
    def Round(self,ind,tag,relabel,thresh,pky,tp,pos):
        #print('ROUND:',ind)
        lab=pky[ind]+tag
        if(lab in tp):
            #print('Already assigned',lab)
            return
        #print('Placing',lab)
        mask=numpy.fabs(pos-pos[ind])<thresh
        m1=tp==relabel
        #print(m1)
        tp[mask*m1]=lab

    #go through peaks, big to small. if it unplace, place pk+a, eg 36i.
    #if any peaks are within thresh in the lists, assign them all to this.
    #if the label is already there, do not place.
    def DoRound(self,tag,relabel,thresh,inty,pky,tp,pos):
        for i in range(len(inty)):  
            self.Round(i,tag,relabel,thresh,pky,tp,pos)

    def DoLeftovers(self,tag,nolabel,pky,tp):
        #check assignments
        left=[]  #is anyone left?
        for pk in numpy.unique(pky):
            if(pk+tag not in tp):
                left.append(pk)
        print('left overs:',left)
        mask=tp==nolabel
        #print('m,asky',mask)
        if(len(tp[mask])==1 and len(left)==1):
            #print('success! we can assign this!')
            tp[mask]=left[0]+tag

    #go through peaks in intensity order, first get
    #unlabelled peaks. if there is no placement of the peak+tag
    #place it. (eg 36(i-1)). Places assignment on most intense
    #peak that is left.
    def DoFinallabel(self,tag,unlabel,pky,tp):
            mask=tp==unlabel #get unassigned peaks left.
            for i in numpy.where(mask==True)[0]:
                 if(pky[i]+tag not in tp): #look for the assignment...
                     tp[i]=pky[i]+tag #assign


    #go over peaks in spec. If there are more peaks than lim (eg 1 for HNCOCA)
    #go over each, find best match to the assigned HNCAs. If the
    #peak number maps to the current peak, save it. otherwise kill it.
    #this is because this peak has been assigned to another slice.
    #we could do the reverse, and place hncoca peaks in other slices.
    #right now we are not doing this....
    def SaveEntriesB(self,ss,tag,spec,lim,inty,pos,tp,stryLine):
        for pk in ss:  #now do the same trick for the HNCOCA
            #self.ShowPeaks('hnca',pk)
            if('hncoca' not in self.peak[pk].keys()): #are there hncoca entries?
                continue
            if(len(self.peak[pk][spec])<=lim): #if there is the right number...
                continue
            for pk3 in self.peak[pk][spec]:  #go over the peaks...
                #go through assignments and figure out who this is.
                vals=[]
                for i in range(len(inty)): #for all the hnca peaks...
                    vals.append(pos[i]-pk3.f3p)  #how close is hncoca peak to hnca
                vals=numpy.array(numpy.fabs(vals)) #get peak that is nearest to hnca
                argy=numpy.argmin(vals)  #get the nearest peak to this
                print(pk3.name,pk3.inty,pk3.f3p,pos[argy],tp[argy])
                if(tp[argy].split(tag)[0]==pk): #if the name corresponds to current guy...
                    continue
                print('hncoca kill:',pk3.name,'from',pk)
                stry='remove,%s,%s,%s # Maps to %s in overlap group %s' % (spec,pk,pk3.name,tp[argy],ss)
                stryLine.append(stry)

    #go over list of peaks in order of intensity
    #are we missing pk+tag?
    #find most intesnse peak for set pk+start (eg 36i)
    #assign the next multiply assigned tag to pk+tag.
    #co means find 'start', then count down. ie start will be more intense than tag
    #co True means start straight away, ie tag will be more intense than start most likey.
    def DoRecover(self,ss,start,tag,inty,pky,tp,co=False):
        for pk in ss:
            #c1= pk+'i' in tp
            c2= pk+tag in tp  #is there an i-1 entry here?
            if(c2==False): #if there is not...
                if(co==False):
                    j=-1  #wait until we've found 'start'
                else:
                    j=1 #start immediately.

                for i in range(len(inty)): #go over the peaks...
                    if(pky[i]==pk): #find peaks that are our guy.
                        #num=len(numpy.where(tp[i]==tp)[0]) #count how many
                        #print('   ',pky[i],name[i],tp[i],inty[i],num)
                        if(j==1):  #if we have already passed the 'i' assignment
                            num=len(numpy.where(tp[i]==tp)[0]) #count how many
                            if(num>1): #this has been multiply assigned already. assign! 
                                tp[i]=pk+tag
                                break
                        if(tp[i]==pk+start): #we have found the 'i' assignment.
                            if(co==False): #if waiting for 'start', signal start now.
                                j=1
                            else:    #otherwise, break.
                                break
            

    #if there is no entry (tag) for a peak.
    #look through options in that slice. assign this to the most intense
    #that has been double counted already. 
    #co flag means start from most intense. co=false means start from i and go down.
    def DoFinalPlace(self,ss,tag,start,inty,tp,co=False):
        tag='(i-1)'
        start='i'
        for pk in ss:
            #c1= pk+'i' in tp
            c2= pk+tag in tp  #is there an i-1 entry here?
            if(c2==False): #if there is not...
                j=1
                for i in range(len(inty)): #go over the peaks...
                    if(pky[i]==pk): #find peaks that are our guy.
                        #num=len(numpy.where(tp[i]==tp)[0]) #count how many
                        #print('   ',pky[i],name[i],tp[i],inty[i],num)
                        if(j==1):  #if we have already passed the 'i' assignment
                            num=len(numpy.where(tp[i]==tp)[0]) #count how many
                            if(num>1): #this has been multiply assigned already. assign! 
                                tp[i]=pk+tag
                                break
                        if(tp[i]==pk+start): #we have found the 'i' assignment.
                            break


    def CheckMixCA(self,ss,speca,specb,thresh,lim,stryLine):  
        print()
        print('checking set...',ss)
        print()

        #get intensities from first spec, and classify as uncagtegories.
        #if any peak is close to second spec, classify with second label. 
        inty,pos,name,pky,tp=self.AssembleSetInfo(ss,speca,specb,thresh,'nocatnocat','(i-1)')
        
        #assignments of peaks
        for i in range(len(inty)):
            print('start:',name[i],inty[i],pos[i],tp[i])    


        #try to place peak labels 'i'. map assignments to all close peaks
        self.DoRound('i','nocatnocat',thresh,inty,pky,tp,pos)

        #anyone unassigned with tag? If there is one leftover
        #and one unassiged, place it.
        self.DoLeftovers('i','nocatnocat',pky,tp)

        #try to place peak labels 'pk+(i-1)'. map assignments to all close peaks labelled (i-1)
        self.DoRound('(i-1)','(i-1)',0.1,inty,pky,tp,pos)
        

        #assignments of peaks
        #for i in range(len(inty)):
        #    print('end:',name[i],inty[i],pos[i],tp[i])    


        #look at remaining unassigned, place label on most intense of who is left.
        self.DoFinallabel('(i-1)','nocatnocat',pky,tp)
        #last step. If there is a nocatnocat, and there is no i-1, place i-1.
        
        #if we are missing an tag (i-1) assignment for a peak, then go through the list,
        #and assign the most intense next mutlpiy assigned peak to be our guy.
        self.DoRecover(ss,'i','(i-1)',inty,pky,tp)
                    

        #assignments of peaks
        for i in range(len(inty)):
            print('end:',name[i],inty[i],pos[i],tp[i])    

        #go through assignments. if the tag contains any labels of the type,
        #save them. generate the correct kill/save labels and save in stryLine
        self.SaveEntries(pky,inty,name,tp,ss,('i','(i-1)'),speca,stryLine)
        
        #last pass: examine the hncoca peaks. how do these line up?
        print()

        self.SaveEntriesB(ss,'(i-1)',specb,lim,inty,pos,tp,stryLine)



 

    def CheckMixCO(self,ss,speca,specb,thresh,lim,stryLine):  
        print()
        print('checking set (CO)...',ss)
        print()


        #assign to LOCAL hnco only.
        #get intensities from first spec, and classify as uncagtegories.
        #if any peak is close to second spec, classify with second label. 
        #if local flag is on, do this only in current slice (don't include others)
        inty,pos,name,pky,tp=self.AssembleSetInfo(ss,speca,specb,thresh,'nocatnocat','(i-1)',local=True)
        
        #assignments of peaks
        for i in range(len(inty)):
            print('start:',name[i],inty[i],pos[i],tp[i])    


        #try to place peak labels 'i'. map assignments to all close peaks
        #take unclassified peaks, starting with most intense, place one per slice.
        self.DoRound('i','nocatnocat',thresh,inty,pky,tp,pos)

        #anyone unassigned with tag? If there is one leftover
        #and one unassiged, place it.
        self.DoLeftovers('i','nocatnocat',pky,tp)

        #try to place peak labels 'pk+(i-1)'. map assignments to all close peaks labelled (i-1)
        #self.DoRound('(i-1)','(i-1)',0.1,inty,pky,tp,pos)
        #self.DoRound('(i-1)','notcatnocat',0.1,inty,pky,tp,pos)
        
        #if there is no entry for (i-1) for a peak.
        #look through options in that slice. assign this to the most intense
        #that has been double counted already. 
        #co flag means start from most intense as i-1 should be more intense than i.
        self.DoRecover(ss,'i','(i-1)',inty,pky,tp,co=True)

        #assignments of peaks
        for i in range(len(inty)):
            print('end:',name[i],inty[i],pos[i],tp[i])    

        #go through assignments. if the tag contains any labels of the type,
        #save them. generate the correct kill/save labels and save in stryLine
        self.SaveEntries(pky,inty,name,tp,ss,('i','(i-1)'),speca,stryLine)
        
        #last pass: examine the hncoca peaks. how do these line up?
        print()

        #self.SaveEntriesB(ss,'(i-1)',specb,lim,inty,pos,tp,stryLine)




    #align spectra and classify peak lists
    def normSpec(self):
        self.ReadManual()

        #or pk in '102','152':
        #    for spec in 'hnco','hncaco':
        #        self.ShowPeaks(spec,pk)
        #for pk in '118','127':
        #    for spec in 'hnca','hncoca':       
        #        self.ShowPeaks(spec,pk)

        #look at overlap in the hnca.
        #use this to reassign overlap resonances.
        #a deficiency: if two peaks are actually close matches,
        #this function will kill them.
        #if there are more than one hncoca, use this to reassign.
        spec='hnca'
        stryLine=[]
        sets=self.GetSets(spec,verb=False)       
        for ss in sets:
            tig=0
            for pk in ss:
                f3p=self.GetPeaksF3p(pk,spec)  
                if(len(f3p)!=2):
                    tig=1
                    break
            if(tig==0):
                continue
            #get peaks from speca. order them by intensity.
            #if peaks are within thresh of specb, classify as i-1.
            #assign unclassified to the peaks in intensity order.
            #if we have one unclassified, and one leftover, assign it.
            #classify the (i-1)s again in order. Any left overs map to 
            #idenfited peak. Assigned (i)s, and (i-1) work out who is kept
            #and who is killed. Go through specB, if there are more peaks
            #than lim, if the peak maps to current peak, keep it, 
            #kill the rest.
            self.CheckMixCA(ss,'hnca','hncoca',0.1,1,stryLine)

        
        self.ManualRemovePeaks(spec='hnco')
        spec='hnco'

        sets=self.GetSets(spec,verb=True) 
        spec='hncaco'      
        for ss in sets:
            tig=0
            for pk in ss:
                f3p=self.GetPeaksF3p(pk,spec)  
                if(len(f3p)>1):
                    tig=1
                    break
            if(tig==0):
                continue
            #get peaks from speca. order them by intensity.
            #if peaks are within thresh of specb, classify as i-1.
            #assign unclassified to the peaks in intensity order.
            #if we have one unclassified, and one leftover, assign it.
            #classify the (i-1)s again in order. Any left overs map to 
            #idenfited peak. Assigned (i)s, and (i-1) work out who is kept
            #and who is killed. Go through specB, if there are more peaks
            #than lim, if the peak maps to current peak, keep it, 
            #kill the rest.
            self.CheckMixCO(ss,'hncaco','hnco',0.1,1,stryLine)



            #now redistribute peaks
                    
                #inty=self.GetPeaksInty(pk,spec)
                #f3p=self.GetPeaksF3p(pk,spec)  
                #self.ShowPeaks(spec,pk)
                #self.ShowPeaks('hncoca',pk)
            #take most intense peak.
            #place this as i in the relevant frame.
            #should place 158i (56)
            #should place 118i (61.3)
            #then 
        if(len(stryLine)!=0):
            for stry in stryLine:
                print (stry)
            #sys.exit(100)
        #sys.exit(100)
        #now repeat the trick for the HNCO/HNCACO.


        print('-------------------------------------------')
        print('Classifying peak lists')
        
        self.AddManualPeaks()
        
        self.ManualRemovePeaks(spec='hncaco')

        specs=self.spec.keys() #classify
        
        self.normCO()  #ifHNCO/HNCACO are present, classify
        self.normCA()
        self.normCACB()
        if('cbcaconh' in specs): #should be able to treat the same as CACB? check.
            if('hnca' in specs):
                self.normCBCACO()
        if('hncocanh' in specs): #adjust charlie's first go.
            if('hncanh' in specs):
                self.normNHNH()
        
        

        self.normTOCSY()
        self.normHTOCSY()


        
        self.SwapPeakClassifications()
        self.ManualRemovePeaks()
        self.ChangePeakLabel()
        self.SetPeakLabel()
        
        self.ScreenHNCOCAinconsistency()  #if hncoca in peaks then clean it up with this.
        self.CheckCA()     #apply hnca/hncoca rules
        self.CheckTOCSY()  #identify CA from ctocsy, look at hnca/hncoca



        

    def MakeToleranceHist(self):
        self.tolHist={}
        if('hncaco' in self.spec.keys() and 'hnco' in self.spec.keys()):
            self.tolHist['CO_CACO']={}
            self.tolHist['CO_CACO']['raw']=[]
            for pk,specs in self.peak.items():
                if('hncaco' not in specs or 'hnco' not in specs):
                    continue
                if(len(specs['hncaco'])>2 or len(specs['hnco'])!=1): #this is the correct number of HNCACO peaks...
                    continue
                
                co=self.GetPeakTp('hnco',pk,'').f3p    
                cacoPk=self.GetPeakTp('hncaco',pk,'')
                if(cacoPk==False):
                    continue
                caco=cacoPk.f3p
                self.tolHist['CO_CACO']['raw'].append(co-caco)
        if('hnca' in self.spec.keys() and 'hncoca' in self.spec.keys()):
            self.tolHist['CA_COCA']={}
            self.tolHist['CA_COCA']['raw']=[]
            for pk,specs in self.peak.items():
                if('hncoca' not in specs or 'hnca' not in specs):
                    continue
                if(len(specs['hnca'])==0 or len(specs['hncoca'])==0): #this is the correct number of HNCACO peaks...
                    continue
                if(len(specs['hnca'])>2 or len(specs['hncoca'])>1): #this is the correct number of HNCACO peaks...
                    continue
                self.ShowPeaks('hnca',pk)
                caPk=self.GetPeakTp('hnca',pk,'')
                if(caPk==False):
                    continue
                ca=caPk.f3p   
                coca=self.GetPeakTp('hncoca',pk,'').f3p       
                self.tolHist['CA_COCA']['raw'].append(ca-coca)

        self.DoG1histograms('G1edges',self.G1edges)

    def DoG1histograms(self,G1lab,G1edges):
        self.tolHist[G1lab]={}
        self.tolHist[G1lab]['raw']=[]
        #self.tolHist[G1lab+'B']={}
        #self.tolHist[G1lab+'B']['raw']=[]
        for peak,edges in G1edges.items(): #all reciprocated edges
            for pk,val,d in edges:
                #if(d=='f'):
                #    self.tolHist[G1lab+'F']['raw'].append(val)
                #elif(d=='b'):
                self.tolHist[G1lab]['raw'].append(val)

        #self.spec_edges_forwards=spec_edges_forwards
        #self.spec_edges_backwards=spec_edges_backwards
        #self.spec_tols_forwards=spec_tols_forwards
        #self.spec_tols_backwards=spec_tols_backwards
        #self.spec_edges_ref=spec_edges_ref
        
        for i,lab in enumerate(self.spec_edges_ref):
            #print(i,lab)
            self.tolHist[G1lab+lab]={}
            self.tolHist[G1lab+lab]['raw']=[]
            for peak,conns in self.spec_edges_forwards[i].items():
                for (pk,val,d) in conns:  #for all possible connections
                    for pkG,valG,dG in G1edges[peak]: #find ones that are reciprocated (G1edges)
                        if(pk==pkG and dG==d):
                            self.tolHist[G1lab+lab]['raw'].append(val)
            """
            self.tolHist[G1lab+lab+'F']={}
            self.tolHist[G1lab+lab+'F']['raw']=[]
            for peak,conns in self.spec_edges_forwards[i].items():
                for (pk,val,d) in conns:
                    for pkG,valG,dG in G1edges[peak]:
                        if(pk==pkG and dG=='f' and d=='f'):
                            self.tolHist[G1lab+lab+'F']['raw'].append(val)

            self.tolHist[G1lab+lab+'B']={}           
            self.tolHist[G1lab+lab+'B']['raw']=[]
            for peak,conns in self.spec_edges_backwards[i].items():
                for (pk,val,d) in conns:
                    for pkG,valG,dG in G1edges[peak]:
                        if(pk==pkG and dG=='b' and d=='b'):
                            self.tolHist[G1lab+lab+'B']['raw'].append(val)
            """            
         
                
    def WriteTolHist(self):
        for key,dicty in self.tolHist.items():
            hist,edges=numpy.histogram(dicty['raw'])
            edges=(edges[:-1]+edges[1:])*0.5
            self.tolHist[key]['hist']=edges,hist
        #outy=open('outy.out','w')
        #for key,dicty in self.tolHist.items():
        #    X,Y= dicty['hist']
        #    for x,y in zip(X,Y):
        #        outy.write('%s\t%f\t%f\n' % (key,x,y))
        #    outy.write('\n\n')
        #outy.close()
        #sys.exit(100)


    def ScreenHNCOCAinconsistency(self):
        stryLine=[]
        #check for consistency.
        for pk,specs in self.peak.items():
            self.GetShift(pk,tocsy=True)
            
            for key,vals in self.shufty.items(): #for each type...
                if(key==''):
                    continue
                if(len(vals)<=1):  #for each peak...
                    continue

                
                vols=[] #extract shifts
                shift={}
                shiftVal={}
                for val in vals: #for each shift...
                    if(val[1] not in shift.keys()):
                        shift[val[1]]=0
                        shiftVal[val[1]]=[]
                    shift[val[1]]+=1
                    shiftVal[val[1]].append(val[0])
                    vols.append(val[0])
                if(numpy.std(vols)<1): #overall standard devation is fine.
                    continue
                #print(shift)
                #TWO MODES.
                #1. if one peak in HNCA, find a match to HNCACO. If not, abort.
                if('hncoca' not in shift.keys()):
                    continue
                if('hnca' not in shift.keys()):
                    continue

                #first, if there are multiple hnca options, compare the most intense from non-main hnca and hncoca
                #if these match, keep these two, pitch the rest (keep main from hnca also.)
                if(shift['hnca']!=1):
                    inty=[]
                    loc=[]
                    for i,pk3 in enumerate(self.peak[pk]['hnca']):
                        print(pk3.name,pk3.tp,pk3.inty)
                        if(pk3.tp==''):
                            inty.append(numpy.fabs(pk3.inty))
                            loc.append(i)
                    argy1=numpy.argmax(inty)
                    hncaf3=self.peak[pk]['hnca'][loc[argy1]]
                    inty=[]
                    loc=[]
                    for i,pk3 in enumerate(self.peak[pk]['hncoca']):
                        #if(pk3.tp==''):
                        inty.append(numpy.fabs(pk3.inty))
                        loc.append(i)
                    argy2=numpy.argmax(inty)
                    hncocaf3=self.peak[pk]['hncoca'][loc[argy2]] 
                    if(numpy.fabs(hncocaf3.f3p - hncaf3.f3p)>1):
                        continue
                    #save these two, kill the rest.
                    for i,pk3 in enumerate(self.peak[pk]['hnca']):
                        if(i!=argy1 and pk3.tp!='main'):
                            print('need to kill:','hnca',pk3.name,pk3.tp)
                            stry='remove,hnca,%s,%s #  keeping only max intensity guy: %s' % (pk,pk3.name,hncaf3.name)
                            stryLine.append(stry)   
                    for i,pk3 in enumerate(self.peak[pk]['hncoca']):
                        if(i!=argy2):
                            print('need to kill:','hncoca',pk3.name,pk3.tp)
                            stry='remove,hncoca,%s,%s #  keeping only max intensity guy: %s' % (pk,pk3.name,hncocaf3.name)
                            stryLine.append(stry)   
                    continue

                #now go over possible hncoca peaks, and if one is a good match in chemical shift to hnca
                #keep this one, pitch the rest.
                valNew=[]
                for pk3 in self.peak[pk]['hncoca']:
                    valNew.append(numpy.fabs(pk3.f3p-shiftVal['hnca'][0]))
                argy=numpy.argmin(valNew)
                if(valNew[argy]<1):
                    #one peak matched. deleting the remaining HNCOCA peaks
                    for i,pk3 in enumerate(self.peak[pk]['hncoca']):
                        if(i!=argy):
                            print('need to kill:',pk3.name)
                            stry='remove,hncoca,%s,%s #  %s and HNCA match (%.2f,%.2f)' % (pk,pk3.name,self.peak[pk]['hncoca'][argy].name,self.peak[pk]['hncoca'][argy].f3p,shiftVal['hnca'][0])
                            stryLine.append(stry)
                    continue

                #2. If we have cTOCSY, then check for their agreement, adding HNCACO peaks as neccessary.
                if('ctocsy' not in shift.keys()):
                    continue
                if(shift['ctocsy']!=1):
                    continue

                stdNew=numpy.std((shiftVal['hnca'],shiftVal['ctocsy']))
                if(stdNew>1):
                    continue
                for i,pk3 in enumerate(self.peak[pk]['hncoca']):
                    #print('need to kill:',pk3.name)
                    stry='remove,hncoca,%s,%s #  HNCACO peak at (%.2f) as TOCSY and HNCA match (%.2f,%.2f)' % (pk,pk3.name,pk3.f3p,shiftVal['hnca'][0],shiftVal['ctocsy'][0])
                    stryLine.append(stry)

                #print('ERROR: inconsistent shifts: %s %.2f %.2f %s' % (key,numpy.average(vols),numpy.std(vols),vols))
                #print()

        if(len(stryLine)>0):
            for line in stryLine:
                print (line)
            sys.exit(100)


    #determine all possible connectivies between peaks
    def assSpec(self):
        print('---------------------------------------------------')
        print('Working out possible connectivities between spectra')
        specs=self.spec.keys()
        self.Optedges={}
        if('hnco' in specs):
            if('hncaco' in specs):
                self.AssCO()
        if('hnca' in specs):
            if('hncoca' in specs):
                self.AssCA(True)  #we have hnca and hncoca
            else:
                self.AssCA(False)  #we only have an hnca.
        if('hncacb' in specs):
            if('hncocacb' in specs):
                self.AssCACB()

        if('hncanh' in specs):
            if('hncocanh' in specs):
                self.AssNHNH()

        self.graphSpec()  #calculate graphs for magma
        
       

    def CombineCACB(self,pk,specs,specA,specB):
        edgeNew=[]
        if(specA in specs.keys() and specB in specs.keys()):
            for edge in self.Optedges[pk][specA]:
                for odge in self.Optedges[pk][specB]:
                    if(edge[0]==odge[0] and edge[2]==odge[2]):
                        score=numpy.fabs(edge[1])+numpy.fabs(odge[1])
                        edgeNew.append((edge[0],score,edge[2]))
        else:
            if(specA in specs.keys()): #otherwise take what's there.
                edgeNew=self.Optedges[pk][specA]
            else:
                edgeNew=self.Optedges[pk][specB]
        return edgeNew

    #calculate graphs from option lists
    def graphSpec(self): #calculate graphs

        
        self.ReadShiftx2()

        self.assRef={}
        for pk in self.peak.keys():
            resns,probs=self.CompareShiftx2(pk)
            for i,resn in enumerate(resns):
                if(resn not in self.assRef.keys()):
                    self.assRef[resn]=[]
                self.assRef[resn].append((pk,probs[i]))

        if(os.path.exists('dat')==False):
            os.system('mkdir dat')
        outy=open('./dat/assRef.out','w')
        for shiftxresi in self.shiftxresis:
            lab=str(shiftxresi)+self.shiftx2[shiftxresi]['resn']
            outy.write('%s : ' % lab)
            if(lab in self.assRef.keys()):
                vals=[]
                for opt in self.assRef[lab]:
                    vals.append(opt[1])
                argy=numpy.argsort(vals)
                for i in range(len(argy)):
                    arg=argy[len(argy)-i-1]
                    outy.write('%s %.2f ' % (self.assRef[lab][arg][0],self.assRef[lab][arg][1]))
            outy.write('\n')
        outy.close()

        self.G1edges={}
        self.G1edgesFull={}

        spec_edges_forwards = []
        spec_tols_forwards = []

        spec_edges_backwards = []
        spec_tols_backwards = []
        spec_edges_ref=[]
        
        if('hnco' in self.spec.keys() and 'hncaco' in self.spec.keys()):
            COf_edges = self.UpdateG1edges('COf')  #read in i+1 correlations into G1edges (f)
            COb_edges = self.UpdateG1edges('COb')    #read in i+1 correlations into G1edges (b)
            spec_edges_ref.append('CO')
            
            spec_edges_backwards.append(COb_edges) 
            spec_edges_forwards.append(COf_edges) # 1f
            spec_tols_backwards.append(self.tolHNCO)
            spec_tols_forwards.append(self.tolHNCO)

        #if('hnca' in self.spec.keys() and 'hncoca' in self.spec.keys()):
        if('hnca' in self.spec.keys() ):
        
            spec_edges_ref.append('CA')
            
            CAb_edges = self.UpdateG1edges('CAb') 
            CAf_edges = self.UpdateG1edges('CAf')    
            spec_edges_forwards.append(CAf_edges)  # 2f
            spec_edges_backwards.append(CAb_edges)
            spec_tols_forwards.append(self.tolHNCA)
            spec_tols_backwards.append(self.tolHNCA)

        if('hncacb' in self.spec.keys() and 'hncocacb' in self.spec.keys()):

            spec_edges_ref.append('CB')
            
            CBb_edges = self.UpdateG1edges('CBb') # b

            CBf_edges = self.UpdateG1edges('CBf') # 5f


            spec_edges_backwards.append(CBb_edges)
            spec_tols_backwards.append(self.tolHNCACB)
   

            spec_edges_forwards.append(CBf_edges)
            spec_tols_forwards.append(self.tolHNCACB)
          



        if('hncanh' in self.spec.keys() and 'hncocanh' in self.spec.keys()):

            spec_edges_ref.append('NH1')
            spec_edges_ref.append('NH2')
            

            NHf1_edges = self.UpdateG1edges('NHf1') #f

            NHf2_edges = self.UpdateG1edges('NHf2') #f
            
            NHb1_edges = self.UpdateG1edges('NHb1') # b
            NHb2_edges = self.UpdateG1edges('NHb2') #b
  


            spec_edges_backwards.append(NHb1_edges)
            spec_edges_backwards.append(NHb2_edges)

            spec_tols_backwards.append(self.tolHNCANH)
            spec_tols_backwards.append(self.tolHNCANH)
     

            spec_edges_forwards.append(NHf1_edges) # 6f
            spec_edges_forwards.append(NHf2_edges) # 7f
            spec_tols_forwards.append(self.tolHNCANH)
            spec_tols_forwards.append(self.tolHNCANH)
            

        
        self.G1edges = {}
        #print(len(spec_edges_forwards))
        #print(spec_edges_forwards[0]['4'])
        #print(spec_edges_forwards[1]['4'])
        #print(spec_edges_backwards[0]['43'])
        #print(spec_edges_backwards[1]['43'])
        self.spec_edges_forwards=spec_edges_forwards
        self.spec_edges_backwards=spec_edges_backwards
        self.spec_tols_forwards=spec_tols_forwards
        self.spec_tols_backwards=spec_tols_backwards
        self.spec_edges_ref=spec_edges_ref

        noe_forward_list = self.Collate_NOE_list(spec_edges_forwards, spec_tols_forwards, 'f')
        noe_backward_list = self.Collate_NOE_list(spec_edges_backwards, spec_tols_backwards, 'b')
        
        
        

        # Manually added edges are in a text file and can add these to the list
        if(os.path.exists('ManualAddedEdges.txt')):
            outy=open('ManualAddedEdges.txt')
            for line in outy.readlines():
                test=line.split()
                if(len(test)==0):
                    continue
                line = line.split('\n')[0]
                line = line.split(',')
                peak_i = line[0]
                direction = line[1]
                connection_peak = line[2]
                if(direction=='f'):
                    noe_forward_list[peak_i][connection_peak] = [numpy.ones(len(spec_edges_forwards)), 0.0]
                    
                elif(direction=='b'):
                    noe_backward_list[peak_i][connection_peak] = [numpy.ones(len(spec_edges_backwards)), 0.0]
                else:
                    print('Need direction to be forward/backward (f/b) in ManualEdgesAdded.txt')
                    print('Line in ManualAddedEdges.txt = ', line)
                    exit()
            outy.close()
        else:
            print('No edges are being manually added')
        
        # Manually remove edges are in a text file and can add remove these from the list
        if(os.path.exists('ManualRemovedEdges.txt')):
            outy=open('ManualRemovedEdges.txt')
            for line in outy.readlines():
                line = line.split('\n')[0]
                line = line.split(',')
                peak_i = line[0]
                direction = line[1]
                connection_peak = line[2]
                if(direction=='f'):
                    try:
                        del noe_forward_list[peak_i][connection_peak]
                    except:
                        pass
                elif(direction=='b'):
                    try:
                        del noe_backward_list[peak_i][connection_peak]
                    except:
                        pass
                else:
                    print('Need direction to be forward/backward (f/b) in ManualRemovedEdges.txt')
                    print('Line in ManualAddedEdges.txt = ', line)
                    exit()
            outy.close()
        else:
            print('No edges are being manually removed')

    
        # save lists of all global 2D peaks where a forward or backward connection hasn't been found 
        # these are saved so that can be checked manually later
        self.checkForward = []
        self.checkBackward = []

        # populate checkForward list with peaks that don't have a forward connection
        for peakLabel in noe_forward_list.keys():
            forward_connections = noe_forward_list[peakLabel].keys()
            if(len(forward_connections)>0):
                found_connection=1
            else:
                found_connection=0
            if(found_connection==0):
                self.checkForward.append(peakLabel)

        # populate checkBackward list with peaks that don't have a backward connection
        for peakLabel in noe_backward_list.keys():
            backward_connections = noe_backward_list[peakLabel].keys()
            if(len(backward_connections)>0):
                found_connection=1
            else:
                found_connection=0
            if(found_connection==0):
                self.checkBackward.append(peakLabel)
       
        

        self.CreateG1edges_list(noe_forward_list, 'f')
        self.CreateG1edges_list(noe_backward_list, 'b')

        #print(self.G1edges['4'])
        #print(self.G1edges['41'])
        #sys.exit(100)        

        self.MakeToleranceHist()


        '''if('hnco' in self.spec.keys() and 'hncaco' in self.spec.keys()): #if we have both, we have a graph
            self.UpdateG1edges('hnco')
            self.UpdateG1edges('hncaco')
            
            #self.UpdateG1edges('hncoca')  #read in i+1 correlations into G1edges (f)
            #self.UpdateG1edges('hnca')    #read in i+1 correlations into G1edges (b)

        


        if('hnca' in self.spec.keys() and 'hncoca' in self.spec.keys()):
            #we can use this to pair down the hnco lists
            #if forward equals hnco forward
            #if backward equals hnca backward
            #keep as unique peak
            for pk in self.G1edges.keys():
                edgeNew=[]
                for edge in self.G1edges[pk]:
                    #if(pk==specTest):
                    #    print(edge)
                    if(edge[2]=='f'): #if edge is 'f'
                        spec='hnca' #compare to label 'hnca'
                    else:
                        spec='hncoca' #compare to label 'hncoca'
                    keep=0

                    #if(pk==specTest):
                    #    print(pk,specTest)
                    #    print(spec in self.Optedges[pk].keys())
                    #    #print(spec,self.Optedges[pk][spec])
                    if(spec in self.Optedges[pk].keys()):
                        for odge in self.Optedges[pk][spec]:
                            if(edge[0]==odge[0]):
                                keep=1
                                break
                    #lse:
                    #    keep=1
                    #if(pk==specTest):
                    #    print(keep)
                    if(keep==1):
                        #print(odge[1],edge[1])
                        score=numpy.fabs(edge[1])+numpy.fabs(odge[1])
                        edgeNew.append((edge[0],score,edge[2]))  

                    # if doesn't have connection in hnca/hncoca, there won't be a connection in edgeNew (i.e. edge/connection has been removed from G1edges)      
                    

                self.G1edges[pk]=edgeNew #update edges list

        #print('after hnca', specTest,self.G1edges[specTest])
        #sys.exit(100)
        #return
        #self.UpdateG1edges('hncacbA')  #read in i+1 correlations into G1edges (f)
        #self.UpdateG1edges('hncacbB')    #read in i+1 correlations into G1edges (b)
        #self.UpdateG1edges('hncacbC')  #read in i+1 correlations into G1edges (f)
        #self.UpdateG1edges('hncacbD')    #read in i+1 correlations into G1edges (b)



        if(5==7 and 'hncacb' in self.spec.keys() and 'hncocacb' in self.spec.keys()):
            #the back and forward lists of these two can be combined to make a single list
            #that is probably unique for forward and backward
            G1new={}
            for pk,specs in self.Optedges.items():
                if(pk not in G1new.keys()):
                    G1new[pk]={}
                specA='hncacbA' #backwards correlation
                specB='hncacbB'
                G1new[pk]['b']=self.CombineCACB(pk,specs,specA,specB)
                specA='hncacbC' #forward correlation
                specB='hncacbD'
                G1new[pk]['f']=self.CombineCACB(pk,specs,specA,specB)
            #print('cb dict:',specTest,G1new[specTest])

            """
            for pk,vals in G1new.items():
                if(pk not in self.G1edges.keys()):
                    self.G1edges[pk]=[]
                for val in vals:
                    for edge in G1new[pk][val]:
                        self.G1edges[pk].append(edge)


            self.tolMatch=0.25
            """

            #mergeCACB dict with main one
            for pk,edges in self.G1edges.items():
                if(pk in G1new):
                    edgeNew=[]
                    for edge in edges:
                        keep=0
                        if(edge[2] in G1new[pk].keys()): #if there are entries
                            for odge in G1new[pk][edge[2]]:
                                if(edge[0]==odge[0]): #need to be the same type
                                    edgeNew.append((edge[0],edge[1],edge[2]))
                                    keep=1
                        else: #if no entries of this type in G1new...
                            edgeNew.append(edge)
                    self.G1edges[pk]=edgeNew

        #print('after hncacb',specTest,self.G1edges[specTest])


        if('hncanh' in self.spec.keys() and 'hncocanh' in self.spec.keys()):
            #we can use this to pair down the hnco lists
            #if forward equals hnco forward
            #if backward equals hnca backward
            #keep as unique peak
            for pk in self.G1edges.keys():
                edgeNew=[]
                for edge in self.G1edges[pk]:
                    #if(pk==specTest):
                    #    print(edge)
                    if(edge[2]=='f'): #if edge is 'f'
                        spec='hncocanhA' #compare to label 'hnca'
                    else:
                        spec='hncanhA' #compare to label 'hncoca'
                    keep=0

                    #if(pk==specTest):
                    #    print(pk,specTest)
                    #    print(spec in self.Optedges[pk].keys())
                    #    #print(spec,self.Optedges[pk][spec])
                    if(spec in self.Optedges[pk].keys()):
                        for odge in self.Optedges[pk][spec]:
                            if(edge[0]==odge[0]):
                                keep=1
                                break
                    #lse:
                    #    keep=1
                    #if(pk==specTest):
                    #    print(keep)
                    if(keep==1):
                        #print(odge[1],edge[1])
                        score=numpy.fabs(edge[1])+numpy.fabs(odge[1])
                        edgeNew.append((edge[0],score,edge[2]))  

                    # if doesn't have connection in hnca/hncoca, there won't be a connection in edgeNew (i.e. edge/connection has been removed from G1edges)      
                    

                self.G1edges[pk]=edgeNew #update edges list'''

    def CreateG1edges_list(self, noe_list, direction):

        
        for peak in noe_list.keys():
            if(peak not in self.G1edges.keys()):
                self.G1edges[peak] = []
            for node2, scores in noe_list[peak].items():
                self.G1edges[peak].append((node2, scores[1], direction))
        #print('CURRENMT')
        #for key,vals in self.G1edges.items():
        #    print(key,vals)
        #sys.exit(100)
            

            
        
       
                
                

    def Collate_NOE_list(self, specs, tols, direction):
        noe_list = {}
        for peak in self.Optedges.keys(): #go through all possible  edge options.
            if(peak not in noe_list.keys()):
                noe_list[peak] = {}
            
            #make sure correlation is within tolerance.
            for i,spec in enumerate(specs): #for connections in each spectrum...
                if(peak in spec.keys()):    #for each peak in the spectrum...
                    for j, edge in enumerate(spec[peak]):  #for all of the edges 
                        if(numpy.fabs(edge[1])<=tols[i]):   #if the score is within the tolerance...
                            entry = 1
                            if(edge[2]==direction): 
                                #if(peak=='4'):
                                #    print ('aaaa',peak,edge,tols[i])
                                    
                                self.edit_array(noe_list,edge,peak,specs, i, entry)
                            
                        elif(numpy.fabs(edge[1])>tols[i]):
                            entry = -1
                            if(edge[2]==direction): 
                                self.edit_array(noe_list,edge,peak,specs, i, entry)

            #if(peak=='4'):
            #    print('fg',noe_list['4'])
            #print('A')
            #try:
            #    print(noe_list['72'])
            #except:
            #    pass
            self.prune_edge_list(noe_list)
            self.prune_edge_list(noe_list)
            #print('B')
            #try:
            #    print(noe_list['72'])
            #except:
            #    pass
            

        #print('ff',noe_list['4'])
        return noe_list

    def prune_edge_list(self, edge_list):
        self.BESTONLY=False #take only the best match numerically
        for node1, node2list in edge_list.items():
            if(len(node2list.keys())>0):
                remove_nodes = []
                
                scrTot=[]
                nodes2rem=[]
                for node2, scores in node2list.items():
                    #find max number of matching cross peaks
                    #print (scores)
                    scrTot.append(numpy.sum(scores[0]))
                    #print(node2,scores)
                    nodes2rem.append(node2)
                #print ('scrtot',scrTot)
                
                maxNo=numpy.max(scrTot) #max number of hits
                
                scrTot=numpy.array(scrTot)
                nodes2rem=numpy.array(nodes2rem)
                mask=(scrTot!=maxNo)
                remove_nodes=nodes2rem[mask]
                #print(remove_nodes)
                for val in remove_nodes:
                    del edge_list[node1][val]


            #for node2, scores in node2list.items():
            #    if(-1 in scores[0]): 
            #        remove_nodes.append(node2)
            #    else:
            #        number_of_zeros = 0
            #        for score in scores[0]:
            #            if(score == 0):
            #                number_of_zeros+=1
            #        if(number_of_zeros>2):
            #            remove_nodes.append(node2)
            #for val in remove_nodes:
            #    del edge_list[node1][val]
            
            #section to only select the edges that have an error less than self.tolMatch

            #if(len(node2list.keys())>0):
                if(self.BESTONLY):
                    weight_list = []
                    connections = []
                    for node2,scores in node2list.items():
                        connections.append(node2)
                        weight_list.append(scores[1])
                    

                    min_weight_index = numpy.argmin(weight_list)
                    best_connection = connections[min_weight_index]
                    remove_connections = []
                    for connection in connections:
                        if(connection!=best_connection):
                            remove_connections.append(connection)
                    
                    for bad_connection in remove_connections:
                        del edge_list[node1][bad_connection]
                



    def edit_array(self, noe_list, edge,peak,specs, i, entry):         
        # edit array [[0,0,...], dist], sets the value in list to 1 if found for that spectrum index
        # if 1 is already in the list, it will then combine the distances and edit the distance to the sum
        #edge = peak,score,direction
        #if(peak=='4'):
        #    print ('b',edge)
        if(edge[0] not in noe_list[peak].keys()):
            noe_list[peak][edge[0]] = [numpy.zeros(len(specs)), numpy.fabs(edge[1])] 
            noe_list[peak][edge[0]][0][i] = entry
        else:
            old_occurance, old_distance = noe_list[peak][edge[0]]
            old_occurance[i] = entry
            new_occurance = old_occurance
            new_distance = numpy.fabs(old_distance) + numpy.fabs(edge[1])
            noe_list[peak][edge[0]] = [new_occurance, new_distance]

        #if(peak=='4'):
        #    print ('c',noe_list['4'])
       

    def DoF1180(self,specAdj):
        if(self.spec[specAdj].f1180!='y'):
            return
        print('Applying f1180 to',specAdj)
        for pk,specs in self.peak.items():
            if(specAdj in specs):
                for i,pk3 in enumerate(self.peak[pk][specAdj]):
                    if(pk3.inty<0):
                        self.DoAlias(pk3,specAdj)

    def SetCSPs(self):
        self.cspH=Parse('assignParFile','cspH')
        if(self.cspH==False):
            self.cspH=0.1
        self.cspN=Parse('assignParFile','cspN')
        if(self.cspN==False):
            self.cspN=0.5
        self.cspN=float(self.cspN)
        self.cspH=float(self.cspH)

    #hnco should have only 1 peak
    #most intense hncaco is i+1
    def normCO(self):
        if('hnco' not in self.spec.keys() and 'hncoca' not in self.spec.keys()):
            return
        print('Classifying HNCO and HNCACO')
        #if peaks are negative, alias them.
        #self.DoF1180('hnco')
        #self.DoF1180('hncaco')

        #this is for the ERRORS
        # Checks if there are the right number of peaks in each hnco bore 

        stryLine=[]  #first, if we have more than 1 in an HNCO, report. Also shout if we have zero.
        for key,vals in self.peak.items():
            if('hnco' not in vals.keys() or len(vals['hnco'])==0):

                #if(len(self.peak[key]['hnco'])==0):
                print('shit: no HNCO for ',key)
                print('Implement maximum peak search.')
                X,Y=self.spec['hnco'].GetSlice(key)    
                argy=numpy.argmax(numpy.fabs(Y))
                print(key,X[argy],Y[argy]/self.spec['hnco'].noise)
                #sys.exit(100)
                stry='add,hnco,%s,%.3f,%s    # no HNCO peaks! try slice maximum with S/N %.2f.' % (key,X[argy],Y[argy],Y[argy]/self.spec['hnco'].noise)
                stryLine.append(stry)

                #print('shit!',key)
                #sys.exit(100)

                continue
            if(len(self.peak[key]['hnco'])==1):
                continue
            
            sets=self.GetSets('hnco',verb=False)
            for ss in sets:
                if(key not in ss):
                    continue
                inty=[]
                for i,pk3 in enumerate(self.peak[key]['hnco']):  #get least intense peak.
                    inty.append(pk3.inty)
                argy=numpy.argmin(numpy.fabs(inty))
                pk3=self.peak[key]['hnco'][argy]
                for s in ss:
                    if(s==key): #if matching self name
                        continue
                    for pk4 in self.peak[s]['hnco']: #otherwise mark for kill
                        if(numpy.fabs(pk4.f3p-pk3.f3p)<0.05):
                            #print('  ',pk4.name,pk4.f3p,pk4.f3p-pk3.f3p,pk4.inty)
                            #print('can safely remove',pk3.name)
                            stry='remove,hnco,%s,%s    # too many HNCO peaks, removing because of overlap with %s' % (key,pk3.name,pk4.name)
                            stryLine.append(stry)
        if(len(stryLine)>0):
            for stry in stryLine:
                print(stry)
            sys.exit(100)        
     

        

        for pk,specs in self.peak.items():
            if('hncaco') not in specs:
                continue

            pk3=self.GetPeakMaxInty(pk,'hncaco')    
            pk3.tp='main' #designate most intense 'main'

            if('hnco' not in specs):
                continue

            pk3=specs['hnco'][0] #take the HNCO peak.
            pk4=self.GetPeaksNear(pk,'hncaco',pk3.f3p) #get nearest peak to hnco in hncaco
            if(numpy.fabs(pk4.f3p-pk3.f3p)>self.tolHNCO): #if there is no peak within threshold...
                self.peak[pk]['hncaco'].append(copy.deepcopy(pk3))  #add.
                self.peak[pk]['hncaco'][-1].name+='a'    #adjust name.
                pk4=self.peak[pk]['hncaco'][-1]  #update pk4

            if( len(self.peak[pk]['hncaco'])<=1 ):
                continue
            
            #classify all other than pk4 as main (forward looking, potential matches)
            for pk5 in self.peak[pk]['hncaco']:
                if(pk4.name==pk5.name):
                    pk5.tp=''
                else:
                    pk5.tp='main'
    

        #this seems to throw up errors and is a bad plan
        stryLine=[]
        for pk,specs in self.peak.items():
            if('hncaco' not in specs or 'hnco' not in specs):
                continue
            if(len(specs['hncaco'])<=2): #this is the correct number of HNCACO peaks...
                continue
            #what to do if we hvae too many.
            #first, get hnco peak.
            pk3=self.peak[pk]['hnco'][0] #
            #now, is 'main' the closest to the hnco?
            vals=[]
            for pk4 in self.peak[pk]['hncaco']:
                vals.append(pk3.f3p-pk4.f3p)
            argy=numpy.argmin(numpy.fabs(vals))
            if(self.peak[pk]['hncaco'][argy].tp!='main'):
                continue
            #reclassify this guy as NOT MAIN
            self.peak[pk]['hncaco'][argy].tp=''
            #go over the others and mark the least intense for death.
            inty=[]
            for pk4 in self.peak[pk]['hncaco']:
                inty.append(pk4.inty)
            inty=numpy.array(inty)
            argys=numpy.argsort(numpy.fabs(inty))

            stry='remove,hncaco,%s,%s    # too many HNCACO peaks, removing least intense bad match to HNCO ' % (pk,self.peak[pk]['hncoca'][argys[0]].name,)
            stryLine.append(stry)
            
            #now reclassify


            #print(inty[argys])
                 #if most 
            #print('PROBLEM')
            #print (pk)

        if(len(stryLine)>0):
            for stry in stryLine:
                print(stry)
            sys.exit(100)       
        
        

            

    def normNHNH(self):
        print('Classifying HNCOCANH and HNCANH')

        # Checks if there are the right number of peaks in each hncocanh bore 
        for key,vals in self.peak.items():
            if('hncocanh' in vals.keys()):
                if(len(self.peak[key]['hncocanh'])>2):
                    print('TOO MANY HNCOCANH PEAKS')
                    print('key = ', key)
                    print('Aborting')
                    #exit()

                spec = 'hncocanh'
                vals=[]
                for i,pk3 in enumerate(self.peak[key][spec]):   # loop through peaks in hncocanh bore
                    vals.append(pk3.inty)

                    ## CB fabs 4/11/19
                    # argy=numpy.argmax(vals)
                    
                argy=numpy.argmax(numpy.fabs(vals))
                self.peak[key][spec][argy].tp='plus'     # call max intensity peak 'plus' (i+1)


        for key,vals in self.peak.items():
            if('hncocanh' in vals.keys()):
                diag_found = False
                if(len(self.peak[key]['hncocanh'])==2):
                    spec = 'hncocanh'
                    diagonal_peakindex = False
                    for i, pk3 in enumerate(self.peak[key][spec]):

                        if(numpy.fabs(pk3.f3p-pk3.f2)<self.tolHNCANH):    # see if there is a weak diagonal peak in the spectrum
                            self.peak[key][spec][i].tp='diag'     # call peak with same Nx and Ny chemical shift diagonal peak 'diag'
                            diag_found= True
                            #print('diagonal found')
                
                    
                    # if diagonal peak is most intense, change less intense peak to plus
                    if(diagonal_peakindex!=False):
                        self.peak[key][spec][int(numpy.fabs(1-diagonal_peakindex))].tp='plus'

                if(diag_found==False):      # if diagonal peak not present, we want to create one

                    peak_diag = copy.deepcopy(self.peak[key][spec][0])
                    #peak_diag.pk += 'd'
                    peak_diag.f3 = peak_diag.f2
                    peak_diag.f3p = peak_diag.f2
                    peak_diag.tp = 'diag'
                    self.peak[key][spec].append(peak_diag)
                    #print('added diagonal')
                           

                    #print('ABORTING')
                    #sys.exit(100)
        
        for key,vals in self.peak.items():
            if('hncanh' in vals.keys()):
                if(len(self.peak[key]['hncanh'])>3):    # check to see if there are the right number of peaks in the hncanh bore
                    print('TOO MANY HNCANH PEAKS')
                    print('key = ', key)
                    # print('Aborting')
                    # exit()
                    

        #if peaks are negative, alias them.
        #THIS MIOGHT CAUSE PROBLEMS: ALIASING TAKE NOW TO MUCH EARLIER.
        #THIS COULD DISRUPT THE ABOVE CODE BLOCK
        #NEED TO ANALYSE IN PLACE
        #self.DoF1180('hncocanh')
        #self.DoF1180('hncanh')
        """
        for pk,specs in self.peak.items():      # Check if peaks are aliased in hncocanh
            if('hncocanh' in specs):
                spec='hncocanh'
                for i,pk3 in enumerate(self.peak[pk][spec]):
                    if(pk3.inty<0):
                        self.DoAlias(pk3,spec)
            if('hncanh' in specs):      # Check if peaks are aliased in hncanh
                spec='hncanh'
                for i,pk3 in enumerate(self.peak[pk][spec]):
                    if(pk3.inty<0):
                        self.DoAlias(pk3,spec)
        """

        for pk,specs in self.peak.items():  
            if('hncanh') in specs:  
                spec='hncanh'
                vals=[]
                for i,pk3 in enumerate(self.peak[pk][spec]):    # loop through hncanh peaks in the bore
                    vals.append(pk3.inty)

                argy=numpy.argmax(numpy.fabs(vals))
                
                self.peak[pk][spec][argy].tp='main'     # call max intensity peak 'main'

                ## CB fabs 4/11/19
                # argy=numpy.argmax(vals)
                tolHNCANH=0.2

                found_diag2 = False
                if(len(self.peak[pk][spec])>=2):
                    for i, pk3 in enumerate(self.peak[pk][spec]):
                        diagonal_peakindex1 = False
                        if(numpy.fabs(pk3.f3p-pk3.f2)<tolHNCANH):     # check to see if weak intensity diagonal peak is found
                            self.peak[pk][spec][i].tp='diag'     # call peak with same Nx and Ny chemical shift diagonal peak 'diag', if not already called main
                            if(i==argy): # if most intense peak is diagonal one
                                self.peak[pk][spec][int(numpy.fabs(argy-1))].tp='main'   # call one of the other peaks 'main'

                            found_diag2= True
                            #print('diagonal found')
                # elif(len(self.peak[pk][spec])<2):
                #     for i, pk3 in enumerate(self.peak[pk][spec]):
                #             self.peak[pk][spec][i].tp='unsure'     
                            

                        


                if(found_diag2==False):      # if weak intensity diagonal peak not found, create one

                    peak_diag1 = copy.deepcopy(self.peak[pk][spec][0])
                    #peak_diag1.pk += 'd'
                    peak_diag1.f3 = peak_diag1.f2
                    peak_diag1.f3p = peak_diag1.f2
                    peak_diag1.tp = 'diag'
                    self.peak[pk][spec].append(peak_diag1)
                    #print('added diagonal')


                
                '''if(len(self.peak[pk][spec])==1): # if only 1 peak in the hncanh bore, add hncocanh to hncanh and add 'a' to bore peak name
                    if('hncocanh' in specs):    
                        if(len(self.peak[pk]['hncocanh'])==1):
                            self.peak[pk][spec].append(copy.deepcopy(self.peak[pk]['hncocanh'][0]))
                            self.peak[pk][spec][-1].name+='a'
                            self.peak[pk][spec][-1].tp='''

        for key,vals in self.peak.items():



            if('hncocanh' in vals.keys() and 'hncanh' in vals.keys()):  # check if plus peak (i+1) in hncanh is at the same position of the peak in the hncocanh, and change if not
                found_plus = -1
                for i,peaki in enumerate(vals['hncocanh']):
                    if(peaki.tp=='plus'):
                        found_plus=i

                
                if(found_plus != -1):


                    pk3=vals['hncocanh'][found_plus] #take plus HNCOCANH peak (i+1).
                    ## CB 4/11/19 - moving vn and vm out of the for loop because otherwise one is reset while the other is set.
                    vn=-1
                    vm=-1
                    for i,pk4 in enumerate(self.peak[key]['hncanh']): # loop through hncanh peaks in the bore

                        if(pk4.tp=='main'):
                            vm=numpy.fabs(pk4.f3p-pk3.f3p)  # difference between hncanh and hncocanh(i+1) max intensity bore peaks
                            imm=i
                        elif(pk4.tp!='diag'):
                            vn=numpy.fabs(pk4.f3p-pk3.f3p)  # else if hncanh bore peak not a diagonal peak, work out the difference between the hncanh bore peak and hncocanh bore peak
                            inn=i
             

                    if(vn>vm and vn!=-1 and vm!=-1):    # if difference between (not 'main' or 'diag', i.e. i-1) hncanh peak and hncocanh(i+1) peak > difference between hncanh 'main' peak and hncocanh(i+1) peak
                        
                        print('SWITCHING HNCANH ASSIGNMENT FOR ',key)
                        self.peak[key]['hncanh'][inn].tp='minus'   
                        self.peak[key]['hncanh'][imm].tp='plus'      
                    elif(vm>vn and vn!=-1 and vm!=-1):
                        print('KEEPING HNCANH ASSIGNMENT FOR ',key)
                        self.peak[key]['hncanh'][inn].tp='plus'   
                        self.peak[key]['hncanh'][imm].tp='minus'    

            
    ####################################################
   
    def SwapPeakClassifications(self):
        #if(os.path.exists('ManualPeaksSwapped.txt')):
        #    outy=open('ManualPeaksSwapped.txt','r')
        if('swap' not in self.manual.keys()):
            return
        for test in self.manual['swap']:
                
                if(len(test)!=4):
                    print('Wrong number of entries. cannot swap.')
                    print(test)
                    print('We need spec,peak,pk1.tp,pk2.tp')
                    sys.exit(100)
                    continue    
          
                pk = test[1]
                spec = test[0]
                initial_label = test[2]
                desired_label = test[3]
                if(spec=='hncanh'):
                    indexes = self.peak[pk][spec]
                    for i in range(len(indexes)):
                        if(self.peak[pk][spec][i].tp==initial_label):
                            self.peak[pk][spec][i].tp = desired_label
                else:
                    indexes = self.peak[pk][spec]
                    i1=-1;i2=-1
                    for i in range(len(indexes)):
                        if(self.peak[pk][spec][i].name==initial_label):
                            i1=i
                            lab1=self.peak[pk][spec][i].tp
                            break
                    for i in range(len(indexes)):
                        if(self.peak[pk][spec][i].name==desired_label):
                            i2=i
                            lab2=self.peak[pk][spec][i].tp
                            break
                    if(i1!=-1 and i2!=-1):
                        print('Swapping:',pk,'in',spec,self.peak[pk][spec][i1].name,self.peak[pk][spec][i1].tp,'and',spec,self.peak[pk][spec][i2].name,self.peak[pk][spec][i2].tp)
                        self.peak[pk][spec][i1].tp=lab2
                        self.peak[pk][spec][i2].tp=lab1
                        
                    else:
                        print('Cannot do swap. Check labels.')
                        print(test)
                        self.ShowPeak(self.peak[pk][spec])
                        sys.exit(100)
       

      

    
    def ManualRemovePeaks(self,spec=False):
        if('remove' not in self.manual.keys()):
            return
        testNew=[]
        for test in self.manual['remove']:
                if(len(test)<3):
                    print('Wrong number of entries. cannot remove.')
                    print(test)
                    print('We need spec,peak,pk1.tp')
                    continue    
       
                #line=line.split('\n')[0].split(',')
                pk=test[1]
                sp=test[0]
                if(spec!=False and spec!=sp):
                    testNew.append(test) #not deleting.
                    continue

                #trying to delete.
                name=test[2]
                #indexes = self.peak[peak][spec]
                #print(self.peak.keys())
                if(pk not in self.peak.keys()):
                    print('cannot find ',pk,'in',self.peak.keys())
                    continue
                if(sp not in self.peak[pk].keys()):
                    print('cannot find ',sp,'in',self.peak[pk].keys())
                    continue

                rem=0
                for i in range(len(self.peak[pk][sp])):
                    print(self.peak[pk][sp][i].name)
                    if(self.peak[pk][sp][i].name==name):
                        print('Removing:',pk,sp,name)
                        del self.peak[pk][sp][i]
                        rem=1
                        break

                if(rem==0):
                    print('Could not remove peak. Check labels.')
                    print(test)
                    self.ShowPeak(self.peak[pk][sp])
                    sys.exit(100)
        self.manual['remove']=testNew  #update list.

                #sys.exit(100)
                #for i in range(len(indexes)-1):
                #    if(self.peak[peak][spec][i].tp==label):
                #        del self.peak[peak][spec][i]
    def ShowPeak(self,arr):  
        for i,pk3 in enumerate(arr):    
            print(pk3.name,pk3.tp,pk3.f3,pk3.f3p,pk3.inty)          

    def ChangePeakLabel(self):
        if('change' not in self.manual.keys()):
            return
        for test in self.manual['change']:
                if(len(test)!=4):
                    print('Wrong number of entries. cannot change.')
                    print(test)
                    print('We need spec,peak,pk1.tp,pk2.tp')
                    continue    
 
                pk=test[1]
                spec=test[0]
                initial_label=test[2]
                final_label=test[3]
                indexes=self.peak[pk][spec]
                ch=0
                for i in range(len(indexes)):
                    if(self.peak[pk][spec][i].tp==initial_label):
                        ch=1
                        self.peak[pk][spec][i].tp=final_label
                        print('Changing:',pk,spec,initial_label,final_label)
                if(ch==0):
                    print('Could not change labels. Check input.')
                    print(test)
                    self.ShowPeak(self.peak[pk][spec])
                    sys.exit(100)

    def SetPeakLabel(self):
        if('set' not in self.manual.keys()):
            return
        for test in self.manual['set']:
                if(len(test)!=3):
                    print('Wrong number of entries. cannot change.')
                    print(test)
                    print('We need spec,peakId,pk2.tp')
                    continue    
 
                name=test[1]
                pk=name.split('_')[0]
                spec=test[0]
                lab=test[2]
                ch=0
                for i,pk2 in enumerate(self.peak[pk][spec]):
                    if(pk2.name==name):
                        ch=1
                        self.peak[pk][spec][i].tp=lab
                        print('Setting:',pk,spec,self.peak[pk][spec][i].tp)    

                if(ch==0):
                    print('Could not change label. Check input.')
                    print(test)
                    self.ShowPeak(self.peak[pk][spec])
                    sys.exit(100)
                #indexes=self.peak[pk][spec]
                #for i in range(len(indexes)):
                #    if(self.peak[pk][spec][i].tp==initial_label):
                #        self.peak[pk][spec][i].tp=final_label
                #        print('Changing:',pk,spec,initial_label,final_label)       


    #####################################################
    def DoAlias(self,pk3,spec):

        smin=self.spec[spec].uc0min
        smax=self.spec[spec].uc0max
        ds=numpy.fabs(self.spec[spec].index0[0]-self.spec[spec].index0[1])
        save=(smin+smax)/2.


        if(pk3.f3>save):
            #alias downwards
            pk3.f3p-=numpy.fabs(smax-smin+ds)
        else:
            #alias upwards
            pk3.f3p+=numpy.fabs(smax-smin+ds)
            if(spec=='hnca' or spec=='hncoca'):
                if(pk3.f3p>70):
                    pk3.f3p-=numpy.fabs(smax-smin+ds)
                    pk3.f3p-=numpy.fabs(smax-smin+ds)
    
    def DoUnAlias_number(self,number,spec):
        
        smin=self.spec[spec].uc0min
        smax=self.spec[spec].uc0max
        folded = number

        if number < smin or number > smax:
            ds=numpy.fabs(self.spec[spec].index0[0]-self.spec[spec].index0[1])
            save=(smin+smax)/2.

            if(number>save):
                #alias downwards
                folded-=numpy.fabs(smax-smin+ds)
            else:
                #alias upwards
                folded+=numpy.fabs(smax-smin+ds)
                if(spec=='hnca' or spec=='hncoca'):
                    if(folded>70):
                        folded-=numpy.fabs(smax-smin+ds)
                        folded-=numpy.fabs(smax-smin+ds)
        return folded


    def normCBCACO(self):


        for pk,specs in self.peak.items():
            if('cbcaconh' in specs and 'hnca' in specs):
                spec='hnca'
                vals = []
                vals_ref=[]
                #print(self.peak[pk])
                for i,pk3 in enumerate(self.peak[pk]['hnca']):
                    #
                    vals_ref.append(pk3.f3p)
                for i,pk3 in enumerate(self.peak[pk]['cbcaconh']):
                    vals.append(pk3.f3)
                for hnca in vals_ref:
                    for i, cbcaconh in enumerate(vals):
                        #print(numpy.abs(cbcaconh-hnca))
                        if numpy.abs(cbcaconh-hnca) < 1.:
                            self.peak[pk]['cbcaconh'][i].tp='high'

                vals_ref = numpy.array(vals_ref)
                vals = numpy.array(vals)

        #exit()

    #most intensit peak is i
    def normCA(self):
        print('Classifying HNCA')

        #self.ShowPeaks('hnca','160')
        #self.ShowPeaks('hncoca','160')

        stryLine=[]
        for pk,specs in self.peak.items():
            if('hnca' not in specs): #classify HNCA main/not main
                continue

            pk3=self.GetPeakMaxInty(pk,'hnca')  #get max intensity from hnca...  
            pk3.tp='main' #designate most intense 'main'

            if('hncoca' not in specs):
                continue

            #1. add the hncoca peak to the HNCA list if only one peak in hncoca and hnca
            if(len(specs['hnca'])==1 and len(specs['hncoca'])==1): #map hncoca onto hnca
                self.peak[pk]['hnca'].append(copy.deepcopy(self.peak[pk]['hncoca'][0]))
                self.peak[pk]['hnca'][-1].name+='a'  
                pk5=self.peak[pk]['hnca'][-1]
                print("ADDED TO HNCA",pk,pk5.name,pk5.f3p)              
                continue
       

            #2. if match between hnca 2nd and hncoca is bad, add 3rd to hnca from the hncoca.
            pk3=self.GetPeakMaxInty(pk,'hncoca')  #get maximum intensity peak in hncoca
            pk4=self.GetPeaksNear(pk,'hnca',pk3.f3p) #get closest HNCA peak to max in HNCACO
            if(pk4.tp=='main' and len(self.peak[pk]['hnca']) > 1):
                self.peak[pk]['hnca'].append(copy.deepcopy(pk3))
                self.peak[pk]['hnca'][-1].name+='a'  
                pk5=self.peak[pk]['hnca'][-1]
                print("ADDED TO HNCA",pk,pk5.name,pk5.f3p)       



            """
            #2. if 'main' in hnca does not match max peak in hncao, switch labels.
            #fails sometimes: remove!
            pk3=self.GetPeakMaxInty(pk,'hncoca')  #get maximum intensity peak in hncoca
            pk4=self.GetPeaksNear(pk,'hnca',pk3.f3p) #get closest HNCA peak to max in HNCACO
            print('a',pk,pk4.name,pk4.tp,pk4.f3p,pk4.inty)
            if(pk4.tp=='main' and len(self.peak[pk]['hnca']) > 1):  #if 'main' is very close to the hncaco peak, then swap labels.
                if(numpy.fabs(pk4.f3p-pk3.f3p)<self.tolHNCA):  #if they are within a ppm...
                    print('b')
                    print(pk,pk4.name,pk4.f3p)
                    pk4.tp==''  #was main, now set to not main.
                    for i,pk5 in enumerate(specs['hnca']): #rename all others to main? bad plan?
                        if(pk5.tp==pk4.tp): #this was in error.
                            pk5.tp='main'
                            print(pk5.name,pk5.f3p,pk5.inty)
                else:
                    #question this. HNCA main is closest to HNCOCA, indicates the OTHER HNCA peak is bad.
                    #adjust second HNCA peak positions with HNCOCA frequncies, and add 'a' as label.
                    if(len(specs['hnca'])!=2): #make sure there are only two peaks in the HNCA
                        continue
                    print('c')
                    for i,pk5 in enumerate(specs['hnca']): #rename all others to main? bad plan?
                        if(pk5.name!=pk4.name):
                            pk5.f3=pk3.f3
                            pk5.f3p=pk3.f3p
                            pk5.name+='a'
            """

            #if(pk=='104'):
            #    sys.exit(100)

        for pk,specs in self.peak.items():
            if('hnca' not in specs): #classify HNCA main/not main
                continue
            if('hncoca' not in self.spec.keys()):
                continue
            #3.if there are 2 peaks in the HNCA, and no entry in HNCOCA, add an ornament in HNCOCA
            if(len(specs['hnca'])==2): #add the hncoca peak to the list if needed
                pk3=self.GetPeakTp('hnca',pk,'') #get main from HNCA
                if('hncoca' not in specs.keys()): #add the peaks to the HNCOCA list, unclassified.
                    self.peak[pk]['hncoca']=[]
                    self.peak[pk]['hncoca'].append(copy.deepcopy(pk3))
                    self.peak[pk]['hncoca'][-1].name+='a'
                    self.peak[pk]['hncoca'][-1].tp=''

        



    def CheckCA(self):
        stryLine=[]
        for pk,specs in self.peak.items():
            if('hnca' not in specs): #classify HNCA main/not main
                continue
            if('hncoca' not in specs):
                continue

            #if there are more than 2 peaks in the HNCA
            #keep main and the one closest to HNCOCA main.
            if(len(self.peak[pk]['hnca'])<=2):
                continue
            if(len(self.peak[pk]['hncoca'])!=1): #if too many 
                continue     

            pk3=self.GetPeakMaxInty(pk,'hncoca') #get the hncoca single peak.

            #get the not-main HNCA peaks
            vals=[]
            loc=[]
            for i,pk4 in enumerate(specs['hnca']):
                if(pk4.tp!='main'):
                    loc.append(i)
                    vals.append(numpy.fabs(pk4.f3p-pk3.f3p))
            argy=numpy.argmin(vals)
            if(vals[argy]<self.tolHNCA): #if match of peaks is within tolerance...
                for i,pk4 in enumerate(specs['hnca']): #we have one good match. kill the others.
                    #print(pk4.name,pk4.tp,pk4.f3p,i,argy,pk4.inty,pk3.name,pk3.f3p,pk3.inty)
                    if(pk4.tp!='main' and i!=loc[argy]):  #kill the others.
                        stry='remove,hnca,%s,%s    # too many HNCA peaks, removing bad matches to the HNCOCA' % (pk,pk4.name,)
                        stryLine.append(stry)
                        continue

            #if we still have a problem, consider the following.
            #this is not a great rule, can lead to errors.
            vals=[]
            loc=[]
            for i,pk4 in enumerate(specs['hnca']):
                if(pk4.tp!='main'):
                    loc.append(i)
                    vals.append(numpy.fabs(pk4.inty))
            argy=numpy.argmax(numpy.fabs(vals))
            locSave=loc[argy]
            testNew=[]
            for i,pk4 in enumerate(specs['hnca']):
                if(pk4.tp!='main' and i!=locSave):
                    stry='remove,%s,%s,%s    # too many HNCA peaks, removing least intense.' % ('hnca',pk,pk4.name,)
                    stryLine.append(stry)

          
            #self.peak[pk]['hnca']=testNew


            
        if(len(stryLine)>0):
            for stry in stryLine:
                print(stry)
            sys.exit(100)
        #sys.exit(100)
            

    #align spectra
    def reference(self):
        print('-------------------------------------------')
        print('Aligning/referencing spectra/peak lists')
        # self.HNCACBmed=self.reffy('hncocacb','max','hncacb','min') #most negative (negmin) should align with most positive in hncocacb
        #self.HNCACBmed=self.reffyCACB('hncocacb','hncacb') ## New referencing for cacb
        
        specSetup='hnco','hncaco','hnca','hncoca','hncacb','ctocsy','hcconh'
        for spec in self.spec.keys():
            if(spec not in specSetup):
                print("Need to setup referencing rules for",spec)
                if(self.spec[spec].ref=='*'):
                    sys.exit(10)
        #self.CBCACONH=self.reffy('cbcaconh','ca','hnca','max')
        self.reffyDefault('hnco')
        self.reffyDefault('hnca')
        self.HNCOmed=self.reffy('hncaco','min','hnco','max') #max in hnco should align with min in hncaco
        self.HNCAmed=self.reffy('hncoca','max','hnca','min') #min in hnca should aign with max in hncoca
        self.HNCAHNCACBmed=self.reffy('hncacb','max','hnca','max') #max in hncacb should match max in hnca
        self.CTOCSYmed=self.reffyCTOCSY('ctocsy','hnca')
        self.HTOCSYmed=self.reffyDefault('hcconh')
        
        self.reffyDefault('cbcaconh')

        print('Setting up f1180:')
        for spec in self.spec.keys():
            self.DoF1180(spec)

        self.Align('1H') #align proton based on median differences in peak lists
        self.Unfold15N()   #try and unfold 15N dimensions, look for common matches plus or minus a few aliases
        self.Align('15N')  #align 15N

            
    def Align(self,nuc):
        #set median H value for each spectrum.
        vals={}
        for pk,specs in self.peak.items():
            for spec,peaks in specs.items():
                if(len(peaks)>0):
                    if(spec not in vals.keys()):
                        vals[spec]={}
                    if(nuc=='1H'):
                        vals[spec][pk]=peaks[0].f1
                    elif(nuc=='15N'):
                        vals[spec][pk]=peaks[0].f2
        ref=vals[self.refSpec]
        del vals[self.refSpec] #get rid of reference from stack

        if(nuc=='1H'):
            self.spec[self.refSpec].Hmed=0
        elif(nuc=='15N'):
            self.spec[self.refSpec].Nmed=0
        for spec,peaks in vals.items(): #for each spectrum that needs Hreferencing
            v=[]
            for pk,f1 in peaks.items(): #for each peak...
                if(pk in ref.keys()):
                    v.append(f1-ref[pk])
            if(nuc=='1H'):
                self.spec[spec].Hmed=numpy.median(v)
            elif(nuc=='15N'):
                self.spec[spec].Nmed=numpy.median(v) 

   

    #take 15N values from peak2D
    #try a few different values of aliasing
    #find common value
    #set all 15N peak.y values to this value.
    def Unfold15N(self):

        for pk in self.spec[self.refSpec].peak2D: #get all reference peak2D places
            spec_resN = []
            current_peaks = []
            for spec,spectra in self.spec.items(): #for all spectra
                
                res=numpy.abs(spectra.index1[1] - spectra.index1[0])
                uc1max=spectra.uc1max  #get Nmax
                uc1min=spectra.uc1min  #get Nmin
                if(uc1max>uc1min):
                    diff=uc1max-uc1min
                else:
                    diff=uc1min-uc1max
                fold=diff + res
                # print(fold_factor)
                multiple = numpy.arange(-2, 3, 1) * fold

                tig=0
                for k, pk2 in enumerate(spectra.peak2D):
                    if pk2.name != pk.name:
                        continue
                    current_peaks.append((multiple+(numpy.ones(5)*(pk2.y-(0)))))  #current location, folded a few times.
                    tig=1
                    break
                if(tig==1): #if we have the peak in this spectrum...
                    spec_resN.append(res)  #save resolution.
            
            if(len(current_peaks)<=1): #if we don't have any values saved...
                continue 

            
            possible_values = current_peaks[0].astype(float) #first row, all possible N chemical shifts...
            if len(possible_values) == 0: #make sure there are possibilities...
                break

            matched=[]
            i = 1
            for row in current_peaks[1:]: #for each other row...
                XX,YY=numpy.meshgrid(possible_values,row) #subtract all possibilities from reference row.
                aa,bb=numpy.where(numpy.fabs(XX-YY)<spec_resN[i]*5)  #find values within a few units of digital resolution...
                if(len(aa)==1): #if we have one match, happy making!
                    for(a,b) in zip(aa,bb):
                        #print('ddd',XX[a,b],YY[a,b])
                        if(i==1):
                            matched.append(YY[a,b])
                        matched.append(XX[a,b])
                i+=1
            print(spec_resN)
            argy=numpy.argmin(spec_resN)
            print(matched)
            print(argy)
            Nval=matched[argy]
            if(numpy.fabs(Nval-numpy.median(matched))>spec_resN[argy]*5):
                print("Problem: referenced value looks problematic")
                print(pk.name,matched)
                sys.exit(100)
                continue

            #update 15N values.
            for spec,spectra in self.spec.items(): #for all spectra...
                for peak in spectra.peak2D:  #for all 2D peaks....
                    if peak.name == pk.name:
                        #if numpy.abs(float(peak.ppmJ-self.ref[spec][1]) - float(possible_values[0])) > 0.001:
                        #    print(peak.name, peak.ppmJ-self.ref[spec][1], possible_values)
                        #if numpy.abs(float(peak.ppmJ) - float(Nval)) > 0.001: #if expected ppm does not match
                        #    print(spec,peak.name, peak.ppmJ, possible_values,Nval)
                        peak.y = Nval  #set unfolded value



    #default referencing: subtract specified number
    def reffyDefault(self,specAdj):
        if(specAdj not in self.spec.keys()):
            return 0
        tol=self.spec[specAdj].ref
        self.AlignPeaks(specAdj,tol)
        return tol
        

    def reffyCTOCSY(self,specAdj,specRef):
        if(specAdj not in self.spec.keys()):
            return 0
        if(specRef not in self.spec.keys()):
            return 0
        #self.HNCACBmed=self.reffy('hncacb','max','hnca','max')

        #Tricky: we don't know which TOCSY peak matched the HNCA.
        #so first, get the CA(i-1) peak from the HNCA.
        #then pick slices with multiple peaks from the TOCSY
        #take all differences and save.
        #then iterate over these lists and find the most common value.
        #set this to CTOCSYmed.
        valCollect=[]
        for pk,specs in self.peak.items():
            if(specAdj in specs and specRef in specs):
                inty=[] #get HNCA main
                for i,pk3 in enumerate(self.peak[pk][specRef]):
                    inty.append(pk3.inty)
                pk3=self.peak[pk][specRef][numpy.argmax(inty)]
                vals=[] #save all TOCSY peaks
                for i,pk2 in enumerate(self.peak[pk][specAdj]):
                    vals.append(numpy.fabs(pk3.f3p-pk2.f3p))
                if(len(vals)<1): #if only 1 tocsy peak, skip.
                    continue
                valCollect.append(numpy.array(vals))
        #print(valCollect)
        vComp=[]  #assumes the FIRST PEAK has the correct difference
        for i,v in enumerate(valCollect[0]): #for each difference encountered
            #print('testing:',v)
            vs=[]
            for j in range(len(valCollect)-1):
                tests=v-valCollect[j+1]
                argNear=numpy.argmin(numpy.fabs(tests))
                vs.append(tests[argNear])
            #print(numpy.average(vs)+v,numpy.std(vs))
            vComp.append((numpy.average(vs)+v,numpy.std(vs))) #get mean and standard deviation
        vComp=numpy.array(vComp)
        #print(vComp)
        argy=numpy.argmin(vComp[:,1])
        tol=vComp[argy,0]
  
        self.AlignPeaks(specAdj,tol)
        return tol
    
    def AlignPeaks(self,specAdj,tol):
        print('Alignment:',specAdj,tol,self.spec[specAdj].ref)
        for pk,specs in self.peak.items(): #adjust peaks
            if(specAdj in specs.keys()):
                for i,pk3 in enumerate(self.peak[pk][specAdj]):
                    self.peak[pk][specAdj][i].f3-=tol
                    self.peak[pk][specAdj][i].f3p-=tol
        self.spec[specAdj].ref=tol

    #return different sorts of max min (abs max and abs min)
    def GetMin(self,ref,vals):
        if(ref=='min'):
            return numpy.argmin(vals)
        if(ref=='max'):
            return numpy.argmax(vals)
        if(ref=='fmin'):
            return numpy.argmin(numpy.fabs(vals))
        if(ref=='fmax'):
            return numpy.argmax(numpy.fabs(vals))

    # COCACB has more than one peak per bore: both peaks should be used during referencing to make sure we get it right.
    def reffyCACB(self, specAdj, specRef):
        
        for pk,specs in self.peak.items():
            if(specAdj in specs):
                vols = []
                dist_COCACB = 0
                for i,pk3 in enumerate(self.peak[pk][specAdj]):
                    vols.append(pk3.f3p)
                if len(vols) == 2:
                    dist_COCACB = vols[0]-vols[1]
            if(specRef in specs):
                vals = []
                intys = []
                dist_CACB = []
                for i,pk3 in enumerate(self.peak[pk][specAdj]):
                    vals.append(pk3.f3p)
                    intys.append(pk3.inty)
                if len(vals) > 1:
                    for x in range(len(vals)) and intys[x] < 0.: ## iterate over peaks and only take the negatives, i.e. cbetas
                        for y in range(len(vals)): ## iterate over all peaks again
                            if x != y and intys[y] > 0.: ## avoid 0 distance peaks and only take positives for distances, i.e. cbeta to calpha distances. 
                                dist_CACB.append(vals[x]-vals[y]) ## append these distances to an array
                #print(dist_COCACB, dist_CACB)
                    
                

    #align adj against ref. peak chosen is either max or min
    #return media tolerance
    def reffy(self,specAdj,optAdj,specRef,optRef):
        if(specAdj not in self.spec.keys()):
            return 0
        if(specRef not in self.spec.keys()):
            return 0
        if(self.spec[specAdj].ref!='*'):
            return self.spec[specAdj].ref
        #compare minimum intensity of CACB with maximum intensity of COCACB
        cnorm=[]    # define new array
        for pk,specs in self.peak.items():     # loop through global 2D peak list (key = peak name), value = list of spectra
            mainy=0
            if(specRef in specs):              # If have reference spectrum
                vols=[]                        # define new array
                for i,pk3 in enumerate(self.peak[pk][specRef]): # loop through list of local peaks
                    vols.append(pk3.inty)       # append peak intensity to vols
                if(len(vols)==0):
                    continue
                argy=self.GetMin(optRef,vols)   # 
                mainy=self.peak[pk][specRef][argy].f3
            if(specAdj in specs and mainy!=0):
                vals=[]
                for i,pk3 in enumerate(self.peak[pk][specAdj]):
                    vals.append(pk3.inty)
                if(len(vals)==0):
                    continue
                argy=self.GetMin(optAdj,vals)
                #argy=numpy.argmax(numpy.fabs(vals))
                cnorm.append(self.peak[pk][specAdj][argy].f3-mainy)
        tol=numpy.median(cnorm)
        
        # flag to remove current referencing problem for hncocacb
        if(specAdj=='hncocacb'):
            tol = 0
        #print('Referencing ',specAdj,' by median: ',tol)
        if(specAdj=='hncocacb'): #NOT SURE WHY I NEED TO DO THIS!
            tol*=-1

        self.AlignPeaks(specAdj,tol) #move peak list, save reference number
        return tol
    
    def normHTOCSY(self):
        if('hcconh' not in self.spec.keys()):
            return
        #no work, boss.
        pass
    def normTOCSY(self):
        if('ctocsy' not in self.spec.keys()):
            return
        #find the CA peak in the list.
        print('Classifying TOCSY')
        stryLine=[]
        for pk,specs in self.peak.items():
            
            if('ctocsy' not in specs):
                continue
            if('hnca' not in specs):
                continue
            go=0
            for i,pk3 in enumerate(self.peak[pk]['hnca']):
                if(pk3.tp!="main"):
                    go=1
                    break
            if(go==0):
                continue
            #find closest peak in tocsy to this

            vals=[]
            for i,pk2 in enumerate(self.peak[pk]['ctocsy']):
                vals.append(numpy.fabs(pk3.f3p-pk2.f3p))
            imin=numpy.argmin(vals)
            if(vals[imin]<1.):
                self.peak[pk]['ctocsy'][imin].tp='CA(i-1)'
                
            
            
    def CheckTOCSY(self):

        print('Checking TOCSY assignment')
        stryLine=[]
        for pk,specs in self.peak.items():
            ca=False
            if('ctocsy' not in specs):
                continue
            if('hnca' not in specs):
                continue
            
            good=0
            
            for i,pk2 in enumerate(self.peak[pk]['ctocsy']):
                #print(pk,pk2.name,pk2.tp,pk2.f3p)
                #print(pk2.tp=='CA(i-1)')
                if(pk2.tp=='CA(i-1)'):
                    good=1
                    break

            #1. if unassigned in the TOCSY.
            #2. if there is one peak in the hnca
            #3. if it agrees well with a cTOCSY peak.
            #4. call the TOCSY (CA(i-1), add an extra peak into the HNCA, classify )
            if(len(self.peak[pk]['hnca'])==1): #if only 1 in HNCA 
                pk3=self.peak[pk]['hnca'][0]  #HNCAO single peak
                pk2=self.GetPeaksNear(pk,'ctocsy',pk3.f3p)  #closest match
                if(numpy.fabs(pk3.f3p-pk2.f3p)<1.): #if we have a match...
                    #we have found our CA(i-1)
                    stry='set,ctocsy,%s,%s   # one HNCA only: matched the TOCSY' % (pk2.name,'CA(i-1)')
                    stryLine.append(stry)
                    #SHOULD PROBABLY TAKE INTENSITY FROM ACTUAL POSITION IN SPECTRUM.
                    stry='add,hnca,%s,%s,0.0   # one HNCA only: bringing in from the TOCSY' % (pk,pk2.f3p)
                    stryLine.append(stry)
                continue
                

            #if(good==1):  #we have a CA assignment. continue.
            #    continue
            #print('no CA assignment:',pk)
            #do some error checking:one test: if we switch the hnca classification, do we get a match?
            #requires hnca, hncoca and ctocsy.
            #take the max intensity HNCA, compare to the TOCSY, if this isn't the case and we match,
            #provide manual strings set this.
            if(len(self.peak[pk]['hnca'])!=2):
                continue

            pk3=self.GetPeakMaxInty(pk,'hnca')    
            if(pk3.tp=='main'): #if the most intense is already main.
                continue

            pk2=self.GetPeaksNear(pk,'ctocsy',pk3.f3p) #get tocsy peak nearest to most intense in HNCA
            if(numpy.fabs(pk3.f3p-pk2.f3p)<1.):
                #we have a match! suggests that we need to change the classification of HNCA.
                #and label the TOCSY.
                #print('SUCCESS')
                pk4=self.GetPeakTp('hnca',pk,'')  #get the 'not main' from hnca

                #pk2.tp='CA(i-1)'
            
                stry='swap,hnca,%s,%s,%s   # taking most intense as main as matches TOCSY' % (pk,pk3.name,pk4.name)
                stryLine.append(stry)
                if(pk2.tp!='CA(i-1)'):
                    stry='set,ctocsy,%s,%s   # matched the TOCSY' % (self.peak[pk]['ctocsy'][imin].name,'CA(i-1)')
                    stryLine.append(stry)

                if('hncoca' not in self.peak[pk].keys()):
                    continue                

                vals=[] #cleanup the HNCOCA
                for i,pk2 in enumerate(self.peak[pk]['hncoca']):
                    #print('hncoca:',pk2.name,pk2.f3,pk4.f3p)
                    test=numpy.fabs(pk4.f3p-pk2.f3p)
                    if(test>1):
                        stry='remove,hncoca,%s,%s   # HNCOCA does not match HNCA' % (pk,pk2.name)
                        stryLine.append(stry)
                        #print(stry)
            
        
        if(len(stryLine)>0):
            for stry in stryLine:
                print(stry)
            #sys.exit(100)
        #sys.exit(100)




    def normCACB(self):
        if('hncacb' not in self.spec.keys() and 'hncocacb' not in self.spec.keys()):
            return 
        print('Classifying CACB:')

        # NOTE: ALIASING NOT DONE
        # file='hncacb2/out/peaksToAlias.txt'
        # listToAlias=[]
        # inny=open(file)
        # for line in inny:
        #     stripped = line.strip('\n')
        #     print(line)
        #     listToAlias.append(stripped)
        for pk,specs in self.peak.items():
            if('hncacb' in specs):
                spec='hncacb'

                vals=[]
                loc=[]
                sgn=numpy.sign(1.) #get sign of positive
                # for i,pk3 in enumerate(self.peak[pk][spec]):
                #     print(pk3.name)
                #     for i,aliasingPeak in enumerate(listToAlias):
                #         if(pk3.name==aliasingPeak):
                #             self.DoAlias(pk3, spec)
                #             print('wahooo')

                #get intensities and locations from hncacb with 'major' sign
                for i,pk3 in enumerate(specs[spec]): 
                    if(sgn==numpy.sign(pk3.inty)):
                        vals.append(pk3.inty)
                        loc.append(i)
                
                if (len(vals)>0): #if we have peaks to analyse...
                    argy=numpy.argsort(-1*numpy.fabs(vals)) #sort big to small         
                    self.peak[pk][spec][loc[argy[0]]].tp='PosMax'  #set most intense to be 'main'
                    if(len(vals)>1): #second most intense, if there, gets minor sign
                        self.peak[pk][spec][loc[argy[1]]].tp='PosMin' 
                #do the same with the sign flipped
                sgn*=-1
                vals=[]
                loc=[]
                for i,pk3 in enumerate(specs[spec]):
                    if(sgn==numpy.sign(pk3.inty)):
                        vals.append(pk3.inty)
                        loc.append(i)
                if (len(vals)>0):
                    argy=numpy.argsort(-1*numpy.fabs(vals)) #sort big to small         
                    sgn=numpy.sign(vals[argy[0]]) #get sign of overall max
                    self.peak[pk][spec][loc[argy[0]]].tp='NegMax'  #set most intense to be 'main'
                    if(len(vals)>1):
                        self.peak[pk][spec][loc[argy[1]]].tp='NegMin'     

            if('hncocacb' in specs):
                spec='hncocacb'
                #will run into trouble if >2 peaks here.
                sgn=numpy.sign(1.0) #take positive sign
                for i,pk3 in enumerate(self.peak[pk][spec]):
                    if(sgn==numpy.sign(pk3.inty)):
                        self.peak[pk][spec][i].tp='Pos'
                    else:
                        self.peak[pk][spec][i].tp='Neg'


            # make sure the label 'negmin' in hncacb is nearest hncocacb
            # if there is already a negmin, then swap if its a posmin or negmax
            #for pk in self.peak.keys():
            if('hncacb' in specs and 'hncocacb' in specs):
                #get hncacb intensities
                vals=[]
                for i,pk2 in enumerate(specs['hncacb']):
                    vals.append(pk2.f3p)
                vals=numpy.array(vals) # ppm of cacb peaks

                for i,pk3 in enumerate(specs['hncocacb']): #get the negative from hncocab
                    if(pk3.tp=='Neg'):
                        break
                argy=numpy.argmin(numpy.fabs(vals-pk3.f3p)) #get closest hncacb to 'neg' in hncocacb
                pk4=self.peak[pk]['hncacb'][argy] #take closest
                #if negmax or posmin and intensity is negative, adjust to negmin.
                if(numpy.fabs(pk4.f3p - pk3.f3p) < self.tolHNCACB): #if difference is within a threshold..
                    if(pk4.tp == 'NegMax'): #this means we have mislabelled hncacb. adjust.
                        print('ADJUSTING CACB LABEL FOR PEAK ',pk4.name,"NEGMIN, was",pk4.name)
                        #print(pk4.tp,pk4.inty)
                        tig=0
                        for i,pk2 in enumerate(specs['hncacb']):
                            if(pk2.tp=='NegMin'):
                                tig=1
                                break
                        if(tig==1):
                            pk2.tp=pk4.tp   #switch the current negmin for pk4's label (NegMax)
                        pk4.tp='NegMin'

    ##############################################
    #routines to get peak information
    def GetPeaksF3p(self,peak,spec,ex=False):
        if(peak not in self.peak.keys()):
            return []
        if(spec not in self.peak[peak].keys()):
            return []
        valR=[]
        for i,pk in enumerate(self.peak[peak][spec]): #take all but the already classified CA
            if(ex==False):
                valR.append(pk.f3p)  
                continue
            if(pk.name[-1]!=ex):
                valR.append(pk.f3p) 

        valR=numpy.array(valR)   
        return valR

    def GetPeaksNear(self,peak,spec,f3):
        if(peak not in self.peak.keys()):
            return []
        if(spec not in self.peak[peak].keys()):
            return []
        valR=[]
        for i,pk in enumerate(self.peak[peak][spec]): #take all but the already classified CA
            valR.append(pk.f3p)  
        valR=numpy.array(valR)
        argy=numpy.argmin(numpy.fabs(valR-f3))
        return self.peak[peak][spec][argy]
        
    def GetPeaksInty(self,peak,spec):
        if(peak not in self.peak.keys()):
            return []
        if(spec not in self.peak[peak].keys()):
            return []
        valR=[]
        for i,pk in enumerate(self.peak[peak][spec]): #take all but the already classified CA
            valR.append(pk.inty)  
        return numpy.array(valR)

    def GetPeakMaxInty(self,peak,spec):    
            vals=self.GetPeaksInty(peak,spec)   #set max intensity CA to main
            argy=numpy.argmax(numpy.fabs(vals))  
            return self.peak[peak][spec][argy]
            

    #show peaks
    def ShowPeaks(self,spec,pk):
        if(pk not in self.peak.keys()):
            return []
        if(spec not in self.peak[pk].keys()):
            return []
        
        
        for pk3 in self.peak[pk][spec]:
            print(spec,pk3.name,pk3.tp,pk3.f3p,pk3.inty)

    #return peak of matching typ
    def GetPeakTp(self,spec,pk,tp):
        if(pk not in self.peak.keys()):
            return False
        if(spec not in self.peak[pk].keys()):
            return False
        
        for pk3 in self.peak[pk][spec]:
            if(pk3.tp==tp):
                return pk3
        return False


    #go through current names, and find one that hasn't been used.
    #ref is the letter used at the end of indexing.
    def FindUnusedName(self,peak,spec,ref):
        go=1
        if(spec not in self.peak[peak].keys()):
            self.peak[peak][spec]=[]

        while(1==1): #find unused name
            tick=0
            for pk3 in self.peak[peak][spec]:
                if(pk3.name[-1]!=ref):
                    continue
                num=pk3.name.split('_')[1][:-1]
                if(num==str(go)):
                    tick=1
            if(tick==0):
                break
            else:
                go+=1
        name=peak+'_'+str(go)+ref     
        return name

    #add a new peak at f3p of type tp to
    #main peak list of peak/spec
    def AddNewPeak(self,peak,spec,f3p,typ,ref):
        #add peak.
        pkl=-1
        for i,n in enumerate(self.spec[spec].peak2D):
            if(n.name==peak):
                pkl=i
                break

        name=self.FindUnusedName(peak,spec,ref)   #get next available name...
        f1=self.spec[spec].peak2D[pkl].x  #get x ppm
        f2=self.spec[spec].peak2D[pkl].y  #get y ppm
        #f3=float(vals[0])                #get z ppm
        f3=f3p
        #f3p=f3                           #initialise aliased z ppm
        inty=self.spec[spec].GetIntensity(peak,f3p)   #get intensity.

        test=name,f1,f2,f3,f3p,inty      #setup line to get peak entry.

        #from .assign_main import peakEntry
        self.peak[peak][spec].append(peakEntry(test))   #add the peak....
        self.peak[peak][spec][-1].tp=typ        #set label
 
    ###############################################

    #take list of chemical shifts, find who is within tolerance and score difference
    def GetPos(self,pk,hnF3,tol):
        #tol = 100
        val=pk.f3p
        mask=numpy.fabs(val-hnF3)<tol #get peaks whose ppm is with tolerance...
        mask2=self.keys[mask]!=pk.pk   #make sure we do not match to ourselves....
        return self.keys[mask][mask2],(hnF3-val)[mask][mask2]  #return the matches.

    def AssCO(self):
        #HNCO is CO(i-1).
        #can get CO(i) by comparing to hncaco major peak (f)
        #can get CO(i-2) by comparing to hncaco main peaks (b)
        self.MatchOptions('COf','hnco','','hncaco','main','f',self.tolHNCO) #look forward
        self.MatchOptions('COb','hncaco','main','hnco','','b',self.tolHNCO) #look backward

    def AssCA(self,hncoca):
        #HNCA main is CA(i).
        #can get CA(i+1) by comparing to hnca minor peak to hnca main peak (f)
        #can get CO(i-1) by comparing to hnca main peak to hncoca peak (or hnca minor peak) (b)
        self.MatchOptions('CAf','hnca','','hnca','main','f',self.tolHNCA) #look forward
        if(hncoca):
            self.MatchOptions('CAb','hnca','main','hncoca','','b',self.tolHNCA) #look backward
        else:
            self.MatchOptions('CAb','hnca','main','hnca','','b',self.tolHNCA) #look backward

    def AssCACB(self):
        #HNCACB is CA(i) and CB(i)
        #can get CA/CB(i-1) by comparing main to minor (b)
        #can get CA/CB(i+1) by comparing minor to major (f)
        # self.MatchOptions('hncacbA','hncacb','PosMax','hncacb','PosMin','b',self.tolHNCACB) #look backward
        #self.MatchOptions('hncacbB','hncacb','NegMax','hncacb','NegMin','b',self.tolHNCACB) #look backward
        # self.MatchOptions('hncacbC','hncacb','PosMin','hncacb','PosMax','f',self.tolHNCACB) #look forward
        #self.MatchOptions('hncacbD','hncacb','NegMin','hncacb','NegMax','f',self.tolHNCACB) #look forward
        self.MatchOptions('CBb','hncacb','NegMax','hncocacb','Neg','b',self.tolHNCACB) #look forward
        self.MatchOptions('CBf','hncocacb','Neg','hncacb','NegMax','f',self.tolHNCACB) #look forward

    def AssNHNH(self):
        self.MatchOptions('NHb1','hncocanh','plus','hncanh','diag','b', self.tolHNCANH)
        self.MatchOptions('NHf1','hncanh','diag','hncocanh','plus','f', self.tolHNCANH)
        self.MatchOptions('NHf2','hncanh','minus','hncocanh','diag','f', self.tolHNCANH)
        self.MatchOptions('NHb2','hncocanh','diag','hncanh','minus','b', self.tolHNCANH)
        #self.MatchOptions('hncanhC','hncanh','diag','hncanh','minus','b', self.tolHNCANH)
        #self.MatchOptions('hncanhD','hncanh','plus','hncanh','diag','b', self.tolHNCANH)



    def MatchOptions(self,lab,refSpec,refTp,matchSpec,matchTp,diry,tol):
        #will define 'lab' in optedges
        #refence list of chemical shifts from refSpec/refTp created.
        #these are then matched to all shifts from matchspec/mathcTp within tolerance.
        #these are used to create an entry in optedges.
        #chemical shift difference and direction are saved

        #first, build shift library indexed to keys
        self.keys=[]
        vals=[]
        for key in self.peak.keys(): #for all peaks in peak list.
            if(refSpec in self.peak[key].keys()):  #if refspec is present.
                for pk2 in self.peak[key][refSpec]:  #for all peaks ..
                    if(pk2.tp==refTp):              #if we match the specified type (refSpec,refTp)
                        self.keys.append(key)       #add the the option list.
                        vals.append(pk2.f3p)
        vals=numpy.array(vals) #chemical shifts of i (CA)
        self.keys=numpy.array(self.keys)
        #now cycle through keys to find matches
        for key in self.peak.keys():                       #for all peaks...
            if(matchSpec in self.peak[key].keys()):        #if the match spectrum is in the list....
                for pk3 in self.peak[key][matchSpec]:      #for all peaks in this list...
                   if(pk3.tp==matchTp):                    #for the combination (matchSpec,matchTp)
                        noe,dist=self.GetPos(pk3,vals,tol)  #for this peak, go through the options, find matches within tolerance.
                        
                        # if (key=='3H-N'):
                        #     print(noe, pk3.f3p)
                        for ii,no in enumerate(noe):  #reveal backward connections. loop over the matches we've just seen...
                            self.UpdateOptedges(lab,pk3,no,dist[ii],diry)  #add entry of type 'lab' in direction 'diry'

        self.SortEdges(lab)

    #sort edge options by chemical shift difference
    def SortEdges(self,spec):
        for key in self.Optedges.keys():
            if(spec in self.Optedges[key]):
                num=[]
                for edge in self.Optedges[key][spec]:
                    num.append(numpy.fabs(edge[1]))  #take the value of the edge...
                argy=numpy.argsort(num)     #sort small to high....
                new=[]  #write out low to high....
                for i in range(len(argy)):  
                    new.append(self.Optedges[key][spec][argy[i]])
                self.Optedges[key][spec]=new  #overwrite sorted values.


    #add edge options into the G1edges graph
    def UpdateG1edges(self,spec):
        G1edges_spec = {}
        for key in self.Optedges.keys():    # for peak in global 2D peak list
            if(spec in self.Optedges[key]):     # if have spectrum present in Optedges for that peak
                if(key not in G1edges_spec.keys()):
                    G1edges_spec[key]=[]
                for edge in self.Optedges[key][spec]:   # for each connection found in OptEdges
                    tig=0
                    for ii,odge in enumerate(G1edges_spec[key]):
                        if(edge[0]==odge[0] and edge[2]==odge[2]):  # if the peak names are identical, and if the connection is in the same direction
                            tig=1
                            break
                    if(tig==0): #edge not here. appending
                        G1edges_spec[key].append(edge)  # add new edge
                    else:
                        scr=numpy.fabs(edge[1])+numpy.fabs(odge[1])
                        #if(key=='4'):
                        #    print('4',edge[0],scr,edge[2] )
                        G1edges_spec[key][ii]=edge[0],scr,edge[2]   # add up distances of two edges (combined error)

        return G1edges_spec
    


    #add options into the connections dictionary
    def UpdateOptedges(self,spec,pk3,no,dist,diry):
        if(pk3.pk not in self.Optedges.keys()):
            self.Optedges[pk3.pk]={}
        if(spec not in self.Optedges[pk3.pk]):
            self.Optedges[pk3.pk][spec]=[]
        self.Optedges[pk3.pk][spec].append( (no,dist,diry ))
        #print(self.Optedges[pk3.pk][spec])

        # (no, dist, diry) = (global_2D_peak_name, how far away peak shift is from exact match peak shift, direction of peak match)




    #verify that connectivities are reciprocal
    def CheckRecip(self,G1edges,key,vals,a,b,strict):
        #print()
        #print('Checking',key,vals,G1edges[key])

        tmp=[]
        scr=[]
        for val in vals: #for each possible edge...
            #print('   ',val[0])
            if(val[2]==a):  #needs to be the correct type...
                if(val[0] in G1edges.keys()):
                    #print('  found',val[0],G1edges[val[0]])
                    for vol in G1edges[val[0]]: #look for reciprocal..
                        #print('  found',val[0],G1edges[val[0]])
                        if(vol[0]==key and vol[2]==b): #if we find it...
                            #print(key,val)
                            #print(val[0],vol,numpy.fabs(vol[1])+numpy.fabs(val[1]))
                            scr.append(numpy.fabs(vol[1])+numpy.fabs(val[1])) #total score.
                            tmp.append(val)  #good good!

        cpy=copy.deepcopy(self.G1edges[key]) #backup old list of connections...

        #print(scr)
        #print(tmp)
        #print(cpy)
        self.G1edges[key]=[]  #blank connections.
        #for tm in tmp:
        #    self.G1edges[key].append(tm)  #append the new options...


        """
        if(strict=='y'):
            if(len(tmp)==1):
                print('keeping:',tmp[0])
                self.G1edges[key].append(tmp[0])
        else:
            if(len(tmp)>0):
                argo=numpy.argmin(scr)
                for i,tm in enumerate(tmp):
                    if(i==argo):
                        print('keeping:',tm)
                        self.G1edges[key].append(tm)

                    elif(numpy.fabs(scr[i])<0.1 ):
                        #elif(numpy.fabs(scr[i]-scr[argo])<0.05 ):
                        self.G1edges[key].append(tm)
                        print('keeping:',tm)
        """
        if(strict=='y'):
            if(len(tmp)>0):
                argo=numpy.argmin(scr)
                for i,tm in enumerate(tmp):
                    if(i==argo):
                        #print(key,'keeping:',tm)
                        self.G1edges[key].append(tm)

                    #elif(numpy.fabs(scr[i])<0.1 ):
                    elif(numpy.fabs(scr[i]-scr[argo])<0.05 ):
                        self.G1edges[key].append(tm)
                        #print(key,'keeping:',tm)
        else:
            if(len(tmp)>0):
                argo=numpy.argmin(scr)
                for i,tm in enumerate(tmp):
                    if(i==argo):
                        #print(key,'keeping:',tm)
                        self.G1edges[key].append(tm)
                    elif(numpy.fabs(scr[i])<self.tolMatch ):
                        #elif(numpy.fabs(scr[i]-scr[argo])<0.05 ):
                        self.G1edges[key].append(tm)
                        #print(key,'keeping:',tm)



        for cp in cpy: #copy in all others that were there previously of the opposite type.
            if(a!=cp[2]):
                self.G1edges[key].append(cp)
                #print(key,'restoring:',cp)
        #return cnt,tmp,scr


    def IsRepeat(self,key):
        fcnt=0 #f count
        bcnt=0 #b count
        rep=0  #repeat count
        #print('looking for repeats:')
        for val in self.G1edges[key]: #for all edges
            if(val[2]=='f'): #if it's forward...
                fcnt+=1
            if(val[2]=='b'): #if it's backward...
                bcnt+=1

            for vol in self.G1edges[key]: #for all other edges..
                if(val!=vol):  #excluding the edge under focus...
                    if(val[0]==vol[0]): #REPEAT
                        if(val[2]=='f'):
                            a=val
                            b=vol
                        else:
                            b=val
                            a=vol
                        rep=1   #increment repeat.
        if(rep==1): #if we have a repeat, and one confident one, remove repeat.
            #print('Repeat found',a,b,fcnt,bcnt)
            if(fcnt+bcnt==3):#if we have only 3 edges to think about
                if(fcnt==2): #if forward  count is 2
                    kl=a
                else:        #if backward count is 2
                    kl=b
                for i,val in enumerate(self.G1edges[key]):
                    if(val==kl): #find matching edge
                        self.G1edges[key].pop(i)
        #print('END:')
        # print(self.G1edges[key])

    def EdgeScreen(self,strict='n'):
        print('-------------------------------------------')
        print('Analysing possible edges. Find edges within')
        print('tolerance, check they are reciprocated.    ')
        #if we see the forward and the backward connection,
        #we have basically sorted this.


        print('GraphStart::')
        for key,vals in self.G1edges.items():

            print (key,vals)
        #print()
        self.G1edgesFull=copy.deepcopy(self.G1edges)

        #specTest='11H-N'
        #print(specTest,self.G1edges[specTest])

        #return
        #sys.exit(100)

        keys=self.G1edges.keys()

        #We do this in three passes.
        #The first takes any reciprocated cross peak pair within combined tolerance tolMatch
        #If 3 are left including a front/back repeat, we get rid of 1 repeat.
        #We then repeat, but take the smallest tol, and  the next if within 0.05
        #We then repeat to clean up the final results and make sure everything is still reciprocated.

        for key in keys:
            vals=self.G1edgesFull[key]
            # print
            # print(key)
            # print(self.G1edges[key])
            self.CheckRecip(self.G1edgesFull,key,vals,'f','b',strict='n')
            self.CheckRecip(self.G1edgesFull,key,vals,'b','f',strict='n')
            self.IsRepeat(key)

        for key in keys:
            vals=self.G1edgesFull[key]
            # print
            # print(key)
            # print(self.G1edges[key])
            self.CheckRecip(self.G1edges,key,vals,'f','b',strict='n')
            self.CheckRecip(self.G1edges,key,vals,'b','f',strict='n')
            self.IsRepeat(key)

        """
        #Take the reciprocated G1edges. 
        #Go over the G1edgesfull
        #purge peaks that have been nailed
        #redo.
        goodList=[]
        self.G1edgesNew={}
        for key,vals in self.G1edges.items():
            if(len(vals)==0):
                continue
            for val in vals:
                pk=val[0]
                dr=val[2]
                if((pk,dr) not in goodList):
                    goodList.append((pk,dr))
        
        print('Assigned:',goodList)
        for key,vals in self.G1edgesFull.items():
            start=self.G1edges[key]         
            for (pk,scr,dir) in vals:
                if((pk,dir) not in goodList):
                    go=1
                    for (pk2,scr2,dir2) in start:
                        if(pk2==pk and dir2==dir):
                            go=0 #already in there
                            break
                    if(go==1):
                        start.append((pk,scr,dir))
            self.G1edgesNew[key]=start

        print('Groph:')
        for key,vals in self.G1edgesFull.items():

            print (key,vals)
        #print()

        for key in keys:
            vals=self.G1edgesNew[key]
            self.CheckRecip(self.G1edgesNew,key,vals,'f','b',strict='n')
            self.CheckRecip(self.G1edgesNew,key,vals,'b','f',strict='n')
            self.IsRepeat(key)

        for key in keys:
            vals=self.G1edgesNew[key]
            self.CheckRecip(self.G1edges,key,vals,'f','b',strict='n')
            self.CheckRecip(self.G1edges,key,vals,'b','f',strict='n')
            self.IsRepeat(key)



        #print('fff')
        print('Graph:')
        for key,vals in self.G1edges.items():

            print (key,vals)
        #print()
        sys.exit(10)
        """
        return

        #for i in range(len(keys)):
        for key in keys:
            #key=keys[argy[i]]
            vals=self.G1edgesFull[key]
            #for key,vals in self.G1edges.items(): #look at all putative connections.
            # print
            # print(key)
            # print(self.G1edges[key])
            #1: look at forwards.
            tp='f' #type of current edge under considation...
            self.CheckRecip(self.G1edges,key,vals,'f','b',strict='y')
            self.CheckRecip(self.G1edges,key,vals,'b','f',strict='y')
            self.IsRepeat(key)


        #for i in range(len(keys)): #do it again to clean up
        for key in keys:
            #key=keys[argy[i]]
            vals=self.G1edgesFull[key]
            #for key,vals in self.G1edges.items(): #look at all putative connections.
            # print
            # print(key)
            # print(self.G1edges[key])
            #1: look at forwards.
            tp='f' #type of current edge under considation...
            self.CheckRecip(self.G1edges,key,vals,'f','b',strict='y')
            self.CheckRecip(self.G1edges,key,vals,'b','f',strict='y')
            self.IsRepeat(key)



    #write inputs for magma
    #also writes an input for mcesCore
    def WriteInit(self,infile=''):

        keys=list(self.G1edges.keys())

        #randomise keys
        #optimisation of magma for backbone
        #goes via a threading algorithm
        vals=range(len(keys))
        from random import shuffle
        shuffle(keys)

        if(infile=='' or self.template==''): #get candidates dictionary
            self.GetConf2()
        else:
            self.GetConf3(infile)



        otty1=open('dat/noe.out','w') #setup data input
        for key in keys:
            for val in self.G1edges[key]:
                otty1.write('%s\t%s\t%s\t%s\n' % (key,val[0],val[1],val[2]))
        otty1.close()

        otty2=open('dat/seq.out','w') #setup sequence input
        for i in range(len(self.G1_nodes)):
            if(self.G1_nodes[i] in self.G1_noes.keys()):
                for val in self.G1_noes[self.G1_nodes[i]]:
                    otty2.write('%s\t%s\t%s\t%s\n' % (self.G1_nodes[i],val[0],val[1],val[2]))
        otty2.close()

        otty3=open('dat/fix.out','w') #setup candidate possibliities
        for key in keys:
            otty3.write('%s : ' % (key))
            for ass in self.candidates[key]:
                otty3.write('%s ' % (ass))
            otty3.write('\n')
        otty3.close()

        """
        #continue to write mcesCore init file
        outy=open('results/mcesCorefull.init','w')
        outy.write('WriteFiles 1\n')
        outy.write('priorIter 1\n') #optimise
        outy.write('maxIter 1000\n')
        outy.write('shortLim 6.0\n')
        outy.write('backbone 1\n')
        outy.write('G1_nodes weight 0\n')
        for key in keys:
            outy.write('%s : ' % key)
            for val in self.G1edges[key]:
                outy.write('%s %.2e %s ' % (val[0],val[1],val[2]))
            outy.write('\n')
        outy.write('G2_nodes weight 0\n')
        for i in range(len(self.G1_nodes)):
            if(self.G1_nodes[i] in self.G1_noes.keys()):
                outy.write('%s : ' % (self.G1_nodes[i]))
                for val in self.G1_noes[self.G1_nodes[i]]:
                    outy.write('%s %s %s ' % (val[0],val[1],val[2]))
                outy.write('\n')
        outy.write('MatchPriorities\n')
        for key in keys:
            outy.write('%s : ' % key)
            for ass in self.candidates[key]:
                outy.write('%s '% (ass))
            outy.write('\n')
        outy.close()
        """

    #read in possible candidates from a file
    def GetConf3(self,infile):
        self.candidates={}
        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            if(len(test)>1):
                tast=line.split(':')
                key=tast[0]
                vals=tast[1].split()
                self.candidates[key]=vals

    #work through residue type possiblities to create candidates dict
    def GetConf2(self):

        self.assRef={}

        self.candidates={}

        
        for key in self.G1edges.keys(): #g1 edges is indexed by peaks
            self.candidates[key]=[]
            resns,probs=self.CompareShiftx2(key)

            

            # print('template = %s' % self.template)
            # print('resns before = %s' % resns)
            if(len(resns)!=0 and self.template!=''):
                # print('have got through if statement 1')
            #if(len(resns)!=0 ):
                for i,resn in enumerate(resns):
                    # print('i = %s' % i)
                    # print('resn = %s' % resn)
                    if(resn[-1]!='P' and int(resn[:-1])!=self.FirstResidue):
                        self.candidates[key].append(resn)
                        # print('updating candidates')
                        # print('self.candidates[key] = %s' % self.candidates[key])
                        if(resn not in self.assRef.keys()):
                            self.assRef[resn]=[]
                        self.assRef[resn].append((key,probs[i] ))
            else:
                for i,G1_node in enumerate(self.G1_nodes): #for each sequence residue
                    #print('G1_node = %s' % G1_node)
                    if(i!=0): #cannot assign to first residue - no NH
                        resn=G1_node[-1] #screen for residue type from sequence name
                        # print(self.shift[key])
                        tig=0
                        if key not in  self.shift.keys():
                            tig = 0
                        else:
                            for koi,vols in self.shift[key]: #if residue type is allowed for i
                                if(koi==resn):
                                    tig=1
                                    # print('tig = 1')
                                    # print('resn1 = %s' %resn)
                                    break

                        if(resn=='P'): #don't add this residue if its a proline
                            tig=0
                        if(tig==1):
                            #if(resn,prob in self.shift[key]): #compare to chemical shift indicies
                            resnim=self.G1_nodes[i-1][-1] #get previous residue in sequence
 
                            # print('resnim = %s' % resnim)
                            tig=0
                            for koi,vils in self.shift2[key]: #go through residues that can be i-1
                                
                                if(koi==resnim):
                                    # print('Have found a match')
                                    tig=1
                                    break
                                # else:
                                #     print('Have not found a match')
                            if(tig==1):  #an assignment option: shift/shift2 both show residue
                                # print('have got through this if statement 2')
                                self.candidates[key].append(G1_node)
                                if(G1_node not in self.assRef.keys()):
                                    self.assRef[G1_node]=[]
                                #print(vols,vils)
                                #print(float(vols)*float(vils))
                                self.assRef[G1_node].append((key,float(vols)*float(vils) ))

            # if(key=='3H-N'):
            #     print('key = %s' % key)
            #     print('self.shift[key] = %s' % self.shift[key])
            #     print('self.shift2[key] = %s' % self.shift2[key])
            #     print('resns = %s' % resns)
            #     print('probs = %s' % probs)
            #     print('candidates = %s' % self.candidates[key])
            #     sys.exit(100)

        #sort assRef by probability
        for key,vals in self.assRef.items():
            srt=[]
            for val in vals:
                srt.append(val[1])
            argy=numpy.flip(numpy.argsort(srt))
            row=[]
            for argy in argy:
                row.append(vals[argy])
            self.assRef[key]=row

        for key,vals in self.candidates.items():
            srt=[]
            for val in vals:
                tig=0
                for koi in self.assRef[val]:
                    if(koi[0]==key):
                        scr=koi[1]
                srt.append(scr)
            argy=numpy.flip(numpy.argsort(srt))
            vals=numpy.array(vals)
            self.candidates[key]=vals[argy]

        #print
        #print(self.assRef)
        #sys.exit(100)
        """
        for key in self.G1edges.keys():
            self.candidates[key]=[]
            for i,G1_node in enumerate(self.G1_nodes): #for each sequence residue

                if(i!=0): #cannot assign to first residue - no NH
                    resn=G1_node[-1] #screen for residue type from sequence name
                    if(resn in self.shift[key]): #compare to chemical shift indicies
                        resnim=self.G1_nodes[i-1][-1]
                        if(resnim in self.shift2[key]):
                            self.candidates[key].append(G1_node)
        """


    #lauch a
    #def launch(self):
    #    import assFrame
    #    assFrame=reload(assFrame)
    #    bool=assFrame.AssFrameMan(self)


    #align peaks in order, taking first occurance of link to be the 'standard'
    def OrderPeaks(self,noes,skip=[]): #set skip to [] to remove orphans, set to False to have everyone
        peaks=[]
        seg=[]
        peakRef=list(self.peak.keys())
        #print(peakRef)
        rem=copy.deepcopy(peakRef)
        while(len(rem)!=0):
            #print ('c')
            curr=rem[0]  
            #print('curr:',curr) 
            peaks,skip=self.AddSeg(peaks,skip,self.GetSegment(curr,'b',noes))
            peaks,skip=self.AddSeg(peaks,skip,self.GetSegment(curr,'f',noes))
            rem=self.GetRemain(peaks,skip,peakRef)
            #print('rem:',rem)  

        #trip the skip list down (can get the seeds from one direction)
        if(skip!=False):
            common=numpy.intersect1d(skip,peaks)   
            skipNew=[]
            for s in skip:
                if(s not in common):
                    skipNew.append(s)
            skip=skipNew
            
            #print(skip)
            #print(peaks)
            #print(len(skip)+len(peaks),len(peakRef))
            
        return numpy.array(peaks),numpy.array(skip)

    def GetRemain(self,peaks,skip,ref):
                if(skip==False):
                    rem=[]
                    for r in ref:
                        if(r not in peaks):
                            rem.append(r)
                    return rem
                #otherwise...
                rem=[]
                for r in ref:
                    if(r not in peaks and r not in skip):
                        rem.append(r)
                return rem

    def GetNext(self,curr,d,noes):
                if(curr not in noes.keys()):
                    return False
                #for noe,(color,val,dir) in noes[curr].items():
                for (noe,val,dir) in noes[curr]:
                        if(dir==d and numpy.fabs(val)<2):
                            return noe
                return False
    def GetSegment(self,curr,d,noes):
                c=curr
                b=[]
                while(c!=False):
                    if(c not in b):
                        b.append(c)
                    else:
                        break
                    c=self.GetNext(c,d,noes) #get next step in direction d.
                if(d=='b'):
                    b=numpy.flip(b)
                return b

    def AddSeg(self,peaks,skip,seg):
                #print('adding seg',peaks,skip,seg)
                if(skip==False):
                    for s in seg:
                        if(s not in peaks):
                            peaks.append(s)
                    return peaks,False
                #otherwise...
                if(len(seg)==1):
                    for s in seg:
                        if(s not in skip):
                            skip.append(s)
                else:
                    for s in seg:
                        if(s not in peaks):
                            peaks.append(s)
              
        
                return peaks,skip

    #calculate gaussian probability of atom,type and measured chemical shift against gross BMRB values
    def DoScr(self,atom,tp,shift):
        if(atom not in self.bmrb.keys()):
            print('atom not in self.bmrb.keys')
            return
        if(tp not in self.bmrb[atom].keys()):
            print('tp not in self.bmrb[atom].keys')
            return

        for val in self.bmrb[atom][tp]: #for all residues of this type....
            lab=self.p3to1[val[0]]  #RESIDUE TYPE
            if(lab not in self.scr.keys()):
                self.scr[lab]={}
            # else:
            #     print('self.scr[lab]')
            #     print(self.scr[lab])
            if(tp not in self.scr[lab].keys()):
                self.scr[lab][tp]={}
            self.scr[lab][tp]=numpy.exp(-(shift-val[1])**2./(2*val[2]**2.))

        # print('self.scr.keys in DoScr')
        # print(self.scr.keys())

    #TOCSY matching function. 
    #nuc can be H or C.
    #shifts=self.bmrbH/C[resn] if doing BMRB
    #shifts=self.shiftx2[resi]['C'/'H'] if doing shiftx
    #if using shiftx, specify cerr.
    #if doing score, and don't just want matches, set SCR to true
    def DoScrTOCSY(self,nuc,peak,bmrb,cerr=False,resn=False,merge=False,IncCA=False):
                
                if(nuc=='C'):
                    spec='ctocsy'
                    STRICT=True #if peak is missing, exclude residue from scr?
                elif(nuc=='H'):
                    spec='hcconh'
                    STRICT=False  #if peak is missing, exclude residue from scr?
                
                #get shifts to use.
                if(peak not in self.peak.keys()):
                    return {}
                if(spec not in self.peak[peak].keys()):
                    return {}
                
                valR=[]
                CAskip=False #will use this later.
                if(nuc=='C'):
                    for i,pk in enumerate(self.peak[peak][spec]): #take all but the already classified CA
                        if(IncCA==True): #if we're taking everything...
                            valR.append(pk.f3p)
                            continue
                        if(pk.tp!='CA(i-1)'):  #otherwise, exclude CA.
                            # print(pk.tp)
                            valR.append(pk.f3p)
                        else:
                            CAskip=True  #CA entry found. don't include CA in list of shifts.
                elif(nuc=='H'):
                    for i,pk in enumerate(self.peak[peak][spec]): #take all but the already classified CA
                        valR.append(pk.f3p)
                if(len(valR)==0): #if we have no matches, give up.
                    return {}
                
                #print('shifts',valR)

                #if we want to loop over all entries...
                #print(resn)
                if(resn==False): #if looping over all residue types...
                    for resn in self.bmrbH.keys(): #these are the amino acids
                        self.DoMatchTOCSY(nuc,valR,resn,bmrb[resn],spec,CAskip=CAskip,SCR=True,cerr=cerr,merge=merge,STRICT=STRICT)
                    return {}
                
                return self.DoMatchTOCSY(nuc,valR,resn,bmrb,spec,CAskip=CAskip,SCR=False,cerr=cerr,merge=merge,STRICT=STRICT) #so strict point to point comparison

    #get shifts from source(bmrb), consolidate if values are closer than specified.
    def ConsolidateShifts(self,spec,nuc,bmrb,cerr=False,merge=False,CAskip=False):
                vol,vop,ti=self.GetShiftArrays(nuc,bmrb,cerr=cerr,CAskip=CAskip)

                #fold peak positions to be within current spectral range.
                for i,vo in enumerate(vol): #for each list from the BMRB...
                    pok=peakEntry(('test','0','0',vo+self.spec[spec].ref,'1')) #create a fake peak entry...
                    self.spec[spec].alias(pok,vo+self.spec[spec].ref,0)  #and alias it to within the range of the spectrum...
                    vol[i]=pok.ppmI-self.spec[spec].ref   #and save.

                #print(merge)
                if(merge==False):
                    return vol,vop,ti
                #if(nuc=='H'):
                #    return vol,vop,ti

                

                #print('   SHIFTS: ',vol,ti)
                #merge values that are close
                #if(nuc=='C'):
                #    close=self.closeC
                #elif(nuc=='H'):
                #    close=self.closeH

                #merge by specified number
                volNew=[]
                vopNew=[]
                tiNew=[]
                vol=numpy.array(vol)
                
                for i in range(len(vol)):
                    diff=numpy.fabs(vol[i]-vol)
                    mask=diff<merge
                    #print(i,vol[i],vol,diff,mask)
                    aa=numpy.where(mask)[0]
                    #print(aa)
                    if(len(aa)==1):
                        volNew.append(vol[i])
                        vopNew.append(vop[i])
                        tiNew.append(ti[i])
                        continue 
                    if(i!=numpy.min(aa)):
                         continue   
                    volNew.append(numpy.average(vol[mask]))
                    vopNew.append(vop[i])

                    tg=ti[i]
                    if(len(ti[i])>2):
                        tg=ti[i][:2]
                    
                    add=tg[-1]
                    for t in ti[mask]:
                        #print(t,tg)
                        
                        if(t[0]==tg[0] and t[1]==tg[1]):
                             pass
                        else:
                            if(t[1] not in add):
                                add+=t[1]
                             #print('SHIT! BAD LABELS!')
                             #print(ti[mask])
                    if(len(add)>1):
                        tg=nuc+add
                        #print(tg)
                        #if(tg!='CBG'):
                        #    sys.exit(100)
                    tiNew.append(tg)
                vol=volNew
                vop=vopNew
                ti=tiNew
                return vol,vop,ti


    def GetShiftArrays(self,nuc,bmrb,cerr=False,CAskip=False):
                typs=bmrb.keys() # these are the Cas, Cbs etc of the given amino acid
                vol=[]
                vop=[]
                ti=[]
                #print(typs)

                for t,vals in bmrb.items():
                        if(t=='C' or t=='H'): #exclude carbonyl and amideH
                            continue

                        try: #unpack different depending on shiftX or bmrb.
                            shift=vals[0]
                            err=vals[1]*1.5
                        except:
                            shift=vals
                            err=cerr
              
                        if(nuc=='C' and shift>110):
                            continue
                        if(nuc=='H' and shift>6):
                            continue
                        
                        if(t=='CA' and CAskip==True): #take all shifts but CA.
                            continue  #skip the CA if we've already identified it.
    
                        vol.append(shift) #exclude carbonyls from TOCSY
                        vop.append(err**2.) # fudge factor for tolerance.
                        ti.append(t) #save type.

                ti=numpy.array(ti)  #save the types as numpy array.
                return vol,vop,ti

    #take expriemntal shifts (valR) and compare to bmrb, both indexed by type.
    #if scoring, populate score dict, and if strict scoring, be harsh.
    #CA is a flag for nuc=c, if there is an indexed CA, then exclude it.
    #can probably do away with this.
    def DoMatchTOCSY(self,nuc,valR,resn,bmrb,spec,CAskip=False,SCR=False,cerr=False,merge=False,STRICT=False):
                
                #get peaks and uncertainties from comparison entry
                if(len(resn)==1): #expand residue to 3 letter to match indexing, if used 1.
                    resn=self.p1to3[resn]
                
                #vol,vop,ti=self.GetShiftArrays(nuc,bmrb,cerr=cerr,CA=CA) #get chemical shift arrays from database.
                
                vol,vop,ti=self.ConsolidateShifts(spec,nuc,bmrb,cerr=cerr,merge=merge,CAskip=CAskip) #get peaks folded into range
                res,cost=self.MatchAtoB(valR,vol,vop,ti,nuc) #compare ValR to Vol, to putatively assign typ.

                if(SCR):  #if scoring, populate scores dictionary.
                    lab=resn
                    if(len(lab)==3):
                        lab=self.p3to1[resn] # residue key/letter
                    if(lab not in self.scr.keys()):
                        self.scr[lab]={}
                    for typ,vals in res.items():
                        if('X' in typ):
                            self.scr[lab][typ]=0
                        else:
                            #res[ti[b]]=valR[a],vol[b],vop[b],prob,a 
                            prob=vals[3] #extract probability...

                return res
                            
                    
        


    def MatchAtoB(self,valR,vol,vop,ti,nuc):
                if(len(valR)==0):
                    return {},-1
                if(len(vol)==0): #make sure residue has at least enough shifts to explain data
                    return {},-1
                
                res={}

                XX,YY=numpy.meshgrid(vol,valR) #create grids of constant row = bmrb and col = expt
                diff=numpy.fabs(XX-YY)**2. #get the absolute value of their differences
                from scipy.optimize import linear_sum_assignment
                ass=linear_sum_assignment(diff)
                
                cost=0
                for a,b in zip(ass[0],ass[1]):
                    #print (a,b)
                    cost+=diff[a,b]**0.5
                    prob=numpy.exp(-1.*(diff[a,b])/(2*vop[b])) #save probability
                    res[ti[b]]=valR[a],vol[b],vop[b],prob,a #save probability
                    #self.peak[peak][spec][vilR[a]].tp=ti[b]+'(i-1)' #classify
                    
                for ii in range(len(valR)): #force unassigned peaks to zero
                    if(ii not in ass[0]):
                        #print('yipyip') #forcing this to zero (crude)
                        res[nuc+'X'+str(ii)]=(valR[ii],) #unassigned and named with an 'X'
                        #if(SCR and STRICT):
                        #    self.scr[lab][nuc+'X'+str(ii)]=0

                return res,cost

    """
    def DoScrTOCSY(self,shifts,nuc='C'):
        #print('----------------------------------------')
        #print('Scoring the %s TOCSY' % nuc)
        if(nuc=='H'):
            valR=[]
            vilR=[]
            #ca=-1
            for i,pk in enumerate(shifts): #take all but the already classified CA
                #if(pk.tp!='CA(i-1)'):
                    # print(pk.tp)
                valR.append(pk.f3p)
                vilR.append(i)
            
            #    else:
            #        ca=pk.f3p
            print('observed shifts:',valR)
            atom='H'  #assuming we've set this to carbon.
            if(atom not in self.bmrb.keys()):
                return

            # print(self.bmrbC.keys())
            #WE NEED TO EXPLAIN ALL PEAKS.

            for resn in self.bmrbH.keys(): #these are the amino acids

                lab=self.p3to1[resn] # residue key/letter
                #print(lab,resn)

                typs=self.bmrbH[resn].keys() # these are the Cas, Cbs etc of the given amino acid
                #print(typs)
                vol=[]
                vop=[]
                ti=[]
                for t in typs:
                    if(t!='H' ): #take all shifts but CA.
                        #if(self.bmrbH[resn][t][0]<120):
                            if(resn=='LYS' and t=='CE'): #HNs on side chain not seen
                                continue
                            vol.append(self.bmrbH[resn][t][0]) #exclude carbonyls from TOCSY
                            vop.append((self.bmrbH[resn][t][1]*1.5)**2.) # fudge factor for tolerance.
                            ti.append(t) #save type.
                ti=numpy.array(ti)
                # print(resn)
                # print(valR) # these are where the peaks are in the tocsy
                # print(vol) # these are where we expect to see the peaks in the tocsy for a given residue

                #reference all BMRB values to within tocsy range (try aliasing)
                #print('pre:',vol)
                for i,vo in enumerate(vol): #for each list from the BMRB...
                    pok=peakEntry(('test','0','0',vo+self.HTOCSYmed,'1')) #create a fake peak entry...
                    self.spec['hcconh'].alias(pok,vo+self.HTOCSYmed,0)  #and alias it to within the range of the spectrum...
                    vol[i]=pok.ppmI-self.HTOCSYmed   #and save.
                #print('vols:',vol)

                #res={}

                if(len(valR)>0 and len(vol)>0): #make sure residue has at least enough shifts to explain data
                    XX,YY=numpy.meshgrid(vol,valR) #create grids of constant row = bmrb and col = expt
                    diff=numpy.fabs(XX-YY)**2. #get the absolute value of their differences
                    from scipy.optimize import linear_sum_assignment
                    ass=linear_sum_assignment(diff)
                    #print('vol',vol)
                    #print('val',valR)
                    #if(len(ass)<len(valR)):
                    #    print('Cannot be this')
                    #print(ass)
                    for a,b in zip(ass[0],ass[1]):
                        #print (a,b)
                        prob=numpy.exp(-1.*(diff[a,b])/(2*vop[b])) #save probability
                        #print('dd',diff[a,b]**0.5,vop[b])
                        #print(ti[b],vol[b],shifts[vilR[a]].f3p,shifts[vilR[a]].tp,prob)
                        if(lab not in self.scr.keys()):
                            self.scr[lab]={}
                        self.scr[lab][ti[b]]=prob #save probability
                    #for i in range(len(valR)): #force unassigned peaks to zero
                    #    #print('looking at ',i,val[i],ass[0])
                    #    if(i not in ass[0]):
                    #        #print('yipyip') #forcing this to zero (crude)
                    #        self.scr[lab]['HX'+str(i)]=0
        

        if(nuc=='C'):
            valR=[]
            vilR=[]
            ca=-1
            for i,pk in enumerate(shifts): #take all but the already classified CA
                if(pk.tp!='CA(i-1)'):
                    # print(pk.tp)
                    valR.append(pk.f3p)
                    vilR.append(i)

                else:
                    ca=pk.f3p
            print('observed shifts:',valR)
            atom='C'  #assuming we've set this to carbon.
            if(atom not in self.bmrb.keys()):
                return

            # print(self.bmrbC.keys())
            #WE NEED TO EXPLAIN ALL PEAKS.


            for resn in self.bmrbC.keys(): #these are the amino acids
                lab=self.p3to1[resn] # residue key/letter
                #print(lab,resn)

                typs=self.bmrbC[resn].keys() # these are the Cas, Cbs etc of the given amino acid
                # print(typs)
                vol=[]
                vop=[]
                ti=[]
                for t in typs:
                    if(t!='CA' and t!='C'): #take all shifts but CA.
                        if(self.bmrbC[resn][t][0]<120):
                            if(resn=='LYS' and t=='CE'): #HNs on side chain not seen
                                continue
                            vol.append(self.bmrbC[resn][t][0]) #exclude carbonyls from TOCSY
                            vop.append((self.bmrbC[resn][t][1]*1.5)**2.) # fudge factor for tolerance.
                            ti.append(t) #save type.
                ti=numpy.array(ti)
                # print(resn)
                # print(valR) # these are where the peaks are in the tocsy
                # print(vol) # these are where we expect to see the peaks in the tocsy for a given residue

                #reference all BMRB values to within tocsy range (try aliasing)
                #print('pre:',vol)
                for i,vo in enumerate(vol): #for each list from the BMRB...
                    pok=peakEntry(('test','0','0',vo+self.CTOCSYmed,'1')) #create a fake peak entry...
                    self.spec['ctocsy'].alias(pok,vo+self.CTOCSYmed,0)  #and alias it to within the range of the spectrum...
                    vol[i]=pok.ppmI-self.CTOCSYmed   #and save.
                #print('vols:',vol)

                #res={}

                if(len(valR)>0 and len(vol)>0): #make sure residue has at least enough shifts to explain data
                    XX,YY=numpy.meshgrid(vol,valR) #create grids of constant row = bmrb and col = expt
                    diff=numpy.fabs(XX-YY)**2. #get the absolute value of their differences
                    from scipy.optimize import linear_sum_assignment
                    ass=linear_sum_assignment(diff)
                    #print('vol',vol)
                    #print('val',valR)
                    #if(len(ass)<len(valR)):
                    #    print('Cannot be this')
                    #print(ass)
                    for a,b in zip(ass[0],ass[1]):
                        #print (a,b)
                        prob=numpy.exp(-1.*(diff[a,b])/(2*vop[b])) #save probability
                        #print('dd',diff[a,b]**0.5,vop[b])
                        #print(ti[b],vol[b],shifts[vil[a]].f3p,shifts[vil[a]].tp,prob)
                        if(lab not in self.scr.keys()):
                            self.scr[lab]={}
                        self.scr[lab][ti[b]]=prob #save probability
                    for i in range(len(valR)): #force unassigned peaks to zero
                        #print('looking at ',i,val[i],ass[0])
                        if(i not in ass[0]):
                            #print('yipyip') #forcing this to zero (crude)
                            self.scr[lab]['CX'+str(i)]=0
    """

                
       
    """
            if(len(valR)>0 and len(vol)>0): #make sure residue has at least enough shifts to explain data
                val=copy.deepcopy(valR)
                vil=copy.deepcopy(vilR)

                XX,YY=numpy.meshgrid(val,vol) #create grids of constant row = bmrb and col = expt
                diff=numpy.fabs(XX-YY) #get the absolute value of their differences
                diffCpy=copy.deepcopy(diff)
                vopCpy=copy.deepcopy(vop)
                tiCpy=copy.deepcopy(ti)

                go=1
                while(go==1):

                    argy=numpy.unravel_index(diff.argmin(), diff.shape) #get minimumum comparison
                    #print('best match:',self.p3to1[resn],ti[argy[0]],vol[argy[0]],shifts[vil[argy[1]]].f3p,shifts[vil[argy[1]]].tp)
                    res[argy[1]]=argy[0] #save assignment

                    if(lab not in self.scr.keys()):
                        self.scr[lab]={}

                    prob=numpy.exp(-1.*(diff[argy])**2./(2*vop[argy[0]]))

                    # if (extraPks > 0):
                    #     prob=prob*(0.8**extraPks)



                    self.scr[lab][ti[argy[0]]]=prob #save probability
                    #print(lab,ti[argy[0]],prob)

                    if(val[argy[1]]>ca and ca!=-1):
                        #print('CATEST')
                        valNew=val[argy[1]]-(self.spec['ctocsy'].alias0) #alias measured shift downwards
                        diffNew=numpy.fabs(valNew-vol[argy[0]])
                        prob2=numpy.exp(-1.*(diffNew)**2./(2*vop[argy[0]])) #save probability
                        #print('previous exp:',val[argy[1]])
                        #print('new aliased: ',valNew)
                        #print('comparison: ',vol[argy[0]])
                        #print(valNew,diffNew,prob2)
                        if(prob2>self.scr[lab][ti[argy[0]]]):
                            self.scr[lab][ti[argy[0]]]=prob2
                    #print(self.scr[lab][ti[argy[0]]])


                    #if testshift>ca, then try probability folding shift back round



                    diff=numpy.delete(diff,argy[0],axis=0) #remove assigned option
                    diff=numpy.delete(diff,argy[1],axis=1) #remove test shift
                    ti=numpy.delete(ti,argy[0],axis=0)     #remove bmrb typ
                    vol=numpy.delete(vol,argy[0],axis=0)   #remove bmrb x0
                    vop=numpy.delete(vop,argy[0],axis=0)   #remove bmrb std

                    val=numpy.delete(val,argy[1],axis=0)
                    vil=numpy.delete(vil,argy[1],axis=0)
                    #print len(ti),len(val)
                    sh=diff.shape
                    if(sh[1]>0): #when we don't have enough columns..
                        pass
                    else:
                        go=0
                    if(sh[0]>0):
                        pass
                    else:
                        go=0

                #find best match for each residue
                #can result in overwriting if multiple assignments hit
                #for ii in range(len(val)):
                #    arg=numpy.argmin(diff[:,ii])
                #    #print(ti[arg],vol[arg],shifts[vil[ii]].f3p,shifts[vil[ii]].tp)
                #    prob=numpy.exp(-1.*(diff[arg,ii])**2./(2*vop[arg]))
                #    if(lab not in self.scr.keys()):
                #        self.scr[lab]={}
                #    self.scr[lab][ti[arg]]=prob #save probability
                
                #from scipy.optimize import linear_sum_assignment
                #ass=linear_sum_assignment(diff)
                #if(len(ass)<len(val)):
                #    print('Cannot be this')
                #print(ass)
                #for a in ass:
                #    print (a)
                #    print(ti[a[0]],vol[a[0]],shifts[vil[a[1]]].f3p,shifts[vil[a[1]]].tp)
                
                print('left over:',val)
                if(len(val)>0): #if we have peaks left over
                    #if peak is at least 5 ppm different from the rest...
                    valR=numpy.array(valR)
                    lab='CX'
                    for vv,v in enumerate(val): #for each remaining peak...
                        doff=numpy.fabs(v-valR)
                        argy=numpy.argsort(doff)
                        print(doff[argy])
                        #if(doff[1]>5): #make sure we are at least 5ppm away from another...
                        if(1==1):
                            #assign it.
                            ii=argy[0]
                            arg=numpy.argmin(diffCpy[:,ii])
                            #print(ti[arg],vol[arg],shifts[vil[ii]].f3p,shifts[vil[ii]].tp)
                            prob=numpy.exp(-1.*(diffCpy[arg,ii])**2./(2*vopCpy[arg]))
                            if(lab not in self.scr.keys()):
                                self.scr[lab]={}
                            if(tiCpy[arg] in self.scr[lab].keys()):
                                print ('AAA')
                                self.scr[lab]['CX'+str(vv+1)]=prob #save probability  
                            else:
                                print('BBBB')
                                self.scr[lab][tiCpy[arg]]=prob #save probability


                #if CB wasnt added, add it with a big negative score.
                if('CB' not in self.scr[lab].keys()):
                    arg=numpy.where(tiCpy=='CB')[0][0]
                    ii=numpy.argmin(diffCpy[arg,:])
                    #print(ti[arg],vol[arg],shifts[vil[ii]].f3p,shifts[vil[ii]].tp)
                    prob=numpy.exp(-1.*(diffCpy[arg,ii])**2./(2*vopCpy[arg]))
                    if(lab not in self.scr.keys()):
                        self.scr[lab]={}
                    self.scr[lab][tiCpy[arg]]=prob #save probability
                continue

                #for ii in range(len(val)):
                #    #argy=numpy.unravel_index(diff.argmin(), diff.shape) #get minimumum comparison
                #    argy=diff[ii].argmin()
                #    print('best match:',val[ii],self.p3to1[resn],ti[argy],vol[argy],shifts[vil[argy[1]]].f3p,shifts[vil[argy[1]]].tp)

 



            else:
                #print 'too many observed shifts. excluding.'
                #if(lab not in self.scr.keys()):
                #    self.scr[lab]={}
                #self.scr[lab]['CG']=0
                pass
    """ 
    def MakeHistogram(self):
        histyRaw={}
        histyRaw['Hhsqc']=[]
        for pk in self.peak.keys():
            for spec in self.peak[pk].keys():
                if(spec not in histyRaw.keys()):
                    histyRaw[spec]=[]
                for i,pk3 in enumerate(self.peak[pk][spec]):
                    #print(pk,spec,i,pk3.f3p)
                    
                    histyRaw[spec].append(pk3.f3p)
                    if(spec=='hncaco'):
                        histyRaw['Hhsqc'].append(pk3.f1)
                
        self.hist={}
        for spec,vals in histyRaw.items():
            gram,edges=numpy.histogram(vals)
            edges=(edges[:-1]+edges[1:])*0.5
            self.hist[spec]=edges,gram
            #print(edges,gram)

        self.MakeBMRBhistogram()


    def MakeBMRBhistogram(self):
        def gaus(x,x0,std):
            return 1/(2*numpy.pi*std**2)**0.5*numpy.exp(-(x-x0)**2./(2*std**2.))

        dat={}
        import pathlib
        bmrbfile=os.path.join(pathlib.Path(__file__).parent.resolve(),'bmrb.txt')
        inny=open(bmrbfile)
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                resn=test[0]
                name=test[1]
                atom=test[2]
                ave=float(test[6])
                std=float(test[7])
                if(resn not in dat.keys()):
                    dat[resn]={}
                if(name not in dat[resn]):
                    dat[resn][name]={}
                dat[resn][name]=ave,std

        def Plot(bins,test,mode=0):
            #binAve=bins[:-1]+(bins[1]+bins[0])/2.
            binAve=bins
            histy=numpy.zeros(len(binAve))
            for resn in dat.keys():
                for name in dat[resn].keys():
                    if(mode==0):
                        if(name==test):
                            #print resn,name
                            histy+=gaus(binAve,dat[resn][name][0],dat[resn][name][1])
                    if(mode==1):
                        if(name[0]==test):
                            #print resn,name
                            histy+=gaus(binAve,dat[resn][name][0],dat[resn][name][1])
            #outy=open('outy.out','a')
            #for i in range(len(binAve)):
            #    outy.write('%e\t%e\n' % (binAve[i],histy[i]))
            #outy.write('\n\n')
            #outy.close()
            return binAve,histy
      
        binNo=200
        #bins=numpy.linspace(4.5,10.,binNo)
        bins=numpy.linspace(0,200.,binNo)
        #bins2D=numpy.linspace(100,150,binNo)     

        self.bmrbHist={}
        self.bmrbHist['CA']=Plot(bins,'CA')
        self.bmrbHist['CB']=Plot(bins,'CB')
        self.bmrbHist['C']=Plot(bins,'C')
        self.bmrbHist['Call']=Plot(bins,'C',mode=1)

        bins=numpy.linspace(-2.5,12.5,binNo)

        self.bmrbHist['HA']=Plot(bins,'HA')
        self.bmrbHist['HB']=Plot(bins,'HB')
        self.bmrbHist['H']=Plot(bins,'H')
        self.bmrbHist['Hall']=Plot(bins,'H',mode=1)

              

    #take spectrum type and peaklabel to work out which atom shift we have
    def GetLab(self,spec,tp):
        #lab=''
        if(spec=='hnco'):
            return 'CO(i-1)'
        if(spec=='hncaco'):
            if(tp=="main"):
                return 'CO(i)'
            else:
                return 'CO(i-1)'
        if(spec=='hnca'):
            if(tp=="main"):
                return 'CA(i)'
            else:
                return 'CA(i-1)'
        if(spec=='hncoca'):
            return 'CA(i-1)'

        if(spec=='hncacb'):
            if(tp=="PosMax"):
                return 'CA(i)'
            elif(tp=="NegMax"):
                return 'CB(i)'
            elif(tp=="PosMin"):
                return 'CA(i-1)'
            elif(tp=="NegMin"):
                return 'CB(i-1)'
        if(spec=='hncocacb'):
            if(tp=="Neg"):
                return 'CB(i-1)'
            elif(tp=='Pos'):
                return 'CA(i-1)'

        if(spec=='cbcaconh'):
            if(tp=="high"):
                return 'CA(i-1)'
            else:
                return 'CB(i-1)'

        if(spec=='hncocanh'):
            if(tp=='plus'):
                return 'N(i+1)'
            if(tp=='diag'):
                return 'N(i)'
        
        if(spec=='hncanh'):
            if(tp=='plus'):
                return 'N(i+1)'
            if(tp=='diag'):
                return 'N(i)'
            if(tp=='minus'):
                return 'N(i-1)'
       
        if(spec=="ctocsy"):
            if(len(tp)>0):
                return tp
        if(spec=="hcconh"):
            if(len(tp)>0):
                return tp

        return ''

        if(spec=="ctocsy"):
            lab=tp

    def CollateErrorsForComboBox(self):
        self.errors_for_combobox = {}
        self.potential_proline_m1 = {}
        self.potential_glycines = {}
        for pk in self.peak.keys():

            if(os.path.exists('CheckedErrorsManually.txt')):
                inny=open('CheckedErrorsManually.txt','r')
                found_peak = False
                for line in inny.readlines():
                    line=line.split('\n')[0]
                    if(line==pk):
                        # Have manually checked the peak and confirmed all is good so don't need to flag it anymore
                        found_peak=True
                inny.close()
                if(found_peak==True):
                    self.errors_for_combobox[pk] = 0
                    self.potential_proline_m1[pk]=0
                    self.potential_glycines[pk] = 0
                    minus1=0
                    Calpha_m1=0
                    if('hncanh' in self.peak[pk].keys() and 'hncocacb' in self.peak[pk].keys()):

                        for i,pk_canh in enumerate(self.peak[pk]['hncanh']):
                            if(pk_canh.tp=='minus'):
                                    minus1+=1
                        for i,pk2 in enumerate(self.peak[pk]['hncocacb']):
                            if(pk2.tp=='Pos'):
                                Calpha_m1=pk2.f3

                        if(minus1==0):
                            if(Calpha_m1!=0 ):
                                if(Calpha_m1>60.0):
                                    self.potential_proline_m1[pk]=1

                        if('hncacb' in self.peak[pk].keys()):
                            nmax=0;nmin=0;pmax=0;pmin=0;
                            for i,pk2 in enumerate(self.peak[pk]['hncacb']):
                                if(pk2.tp=='PosMax'):
                                    pmax+=1
                                if(pk2.tp=='PosMin'):
                                    pmin+=1
                                if(pk2.tp=='NegMax'):
                                    nmax+=1
                                if(pk2.tp=='NegMin'):
                                    nmin+=1

                            if(nmax + nmin < 2):
                                if(pk2.f3 < 47.0):
                                    self.potential_glycines[pk]=1
                                    
                    
                    # Have manually checked the peak and confirmed all is good so don't need to flag it anymore, so continue onto next peak
                    continue
                
                else:

                    self.errors_for_combobox[pk] = 0
                    self.potential_proline_m1[pk]= 0 
                    self.potential_glycines[pk] = 0
                    Calpha_m1=0
                    minus1=0
                    specord='hnco','hncaco','hnca','hncoca','hncacb','hncocacb','hncocanh','hncanh'
                    for spec in specord:
                        if(spec in self.peak[pk].keys()):
                            if(spec=='hnco'):
                                if(len(self.peak[pk][spec])>1):
                                    self.errors_for_combobox[pk]=1
                            if(spec=='hncaco'):
                                if(len(self.peak[pk][spec])>2):
                                    self.errors_for_combobox[pk]=1
                            if(spec=='hnca'):
                                if(len(self.peak[pk][spec])>2):
                                    self.errors_for_combobox[pk]=1
                            if(spec=='hncoca'):
                                if(len(self.peak[pk][spec])>1):
                                    self.errors_for_combobox[pk]=1
                            if(spec=='hncacb'):
                                if(len(self.peak[pk][spec])>4):
                                    self.errors_for_combobox[pk]=1
                                nmax=0;nmin=0;pmax=0;pmin=0;
                                for i,pk2 in enumerate(self.peak[pk][spec]):
                                    if(pk2.tp=='PosMax'):
                                        pmax+=1
                                    if(pk2.tp=='PosMin'):
                                        pmin+=1
                                    if(pk2.tp=='NegMax'):
                                        nmax+=1
                                    if(pk2.tp=='NegMin'):
                                        nmin+=1
                                if(pmin>1 or pmax>1 or nmax>1 or nmin>1):
                                    self.errors_for_combobox[pk]=1
                                if(pmax==1 and pmin==0):
                                    self.errors_for_combobox[pk]=1
                                if(nmax + nmin < 2):
                                    if(pk2.f3 < 47.0):
                                        self.potential_glycines[pk]=1

                            if(spec=='hncocacb'):
                                if(len(self.peak[pk][spec])>2):
                                    self.errors_for_combobox[pk]=1
                                for i,pk2 in enumerate(self.peak[pk][spec]):
                                    if(pk2.tp=='Pos'):
                                        Calpha_m1=pk2.f3
                            
                            if(spec=='hncocanh'):
                                if(len(self.peak[pk][spec])>2):
                                    self.errors_for_combobox[pk]=1
                            
                            if(spec=='hncanh'):
                                if(len(self.peak[pk][spec])>3):
                                    self.errors_for_combobox[pk]=1
                                
                    
                                for i,pk_canh in enumerate(self.peak[pk][spec]):
                                    if(pk_canh.tp=='minus'):
                                        minus1+=1
                    
                    if(minus1==0):
                        if(Calpha_m1!=0 ):
                            if(Calpha_m1>60.0):
                                self.potential_proline_m1[pk]=1

                                
                    
                    self.GetShift(pk)
                    for key,vals in self.shufty.items():
                        if(len(vals)>1):
                            vols=[] #extract shifts
                            for val in vals:
                                vols.append(val[0])
                            if(numpy.std(vols)>1):
                                self.errors_for_combobox[pk]=1

                                diffs=[]
                                for i in range(len(vols)):
                                    diff=0
                                    for j in range(len(vols)):
                                        if(i!=j):
                                            diff+=vols[i]-vols[j]
                                    diffs.append(diff)
                                argy=numpy.argmax(diffs)
                                spoc=vals[argy][1]
                                ii=vals[argy][2]
                                pkErr=self.peak[pk][spoc][ii]
                                self.errors_for_combobox[pk]=1
            
            else:
                self.errors_for_combobox[pk] = 0
                self.potential_proline_m1[pk]= 0
                self.potential_glycines[pk] = 0
                specord='hnco','hncaco','hnca','hncoca','hncacb','hncocacb','hncocanh','hncanh','ctocsy','hcconh'
                Calpha_m1=0
                minus1=0
                for spec in specord:
                    
                    if(spec in self.peak[pk].keys()):
                        if(spec=='hnco'):
                            if(len(self.peak[pk][spec])>1):
                                self.errors_for_combobox[pk]=1
                        if(spec=='hncaco'):
                            if(len(self.peak[pk][spec])>2):
                                self.errors_for_combobox[pk]=1
                        if(spec=='hnca'):
                            if(len(self.peak[pk][spec])>2):
                                self.errors_for_combobox[pk]=1
                        if(spec=='hncoca'):
                            if(len(self.peak[pk][spec])>1):
                                self.errors_for_combobox[pk]=1
                        if(spec=='hncacb'):
                            if(len(self.peak[pk][spec])>4):
                                self.errors_for_combobox[pk]=1
                            nmax=0;nmin=0;pmax=0;pmin=0;
                            for i,pk2 in enumerate(self.peak[pk][spec]):
                                if(pk2.tp=='PosMax'):
                                    pmax+=1
                                if(pk2.tp=='PosMin'):
                                    pmin+=1
                                if(pk2.tp=='NegMax'):
                                    nmax+=1
                                if(pk2.tp=='NegMin'):
                                    nmin+=1
                            if(pmin>1 or pmax>1 or nmax>1 or nmin>1):
                                self.errors_for_combobox[pk]=1
                            if(pmax==1 and pmin==0):
                                self.errors_for_combobox[pk]=1
                            if(nmax + nmin < 2):
                                if(pk2.f3 < 47.0):
                                    self.potential_glycines[pk]=1


                        if(spec=='hncocacb'):
                            if(len(self.peak[pk][spec])>2):
                                self.errors_for_combobox[pk]=1
                            for i,pk2 in enumerate(self.peak[pk][spec]):
                                if(pk2.tp=='Pos'):
                                    Calpha_m1=pk2.f3
                        
                        if(spec=='hncocanh'):
                            if(len(self.peak[pk][spec])>2):
                                self.errors_for_combobox[pk]=1
                        
                        if(spec=='hncanh'):

                            if(len(self.peak[pk][spec])>3):
                                self.errors_for_combobox[pk]=1
                            
                    
                            for i,pk_canh in enumerate(self.peak[pk][spec]):
                                if(pk_canh.tp=='minus'):
                                    minus1+=1
                        
                if(minus1==0):
                    if(Calpha_m1!=0 ):
                        if(Calpha_m1>60.0):
                            self.potential_proline_m1[pk]=1
                
                self.GetShift(pk,tocsy=False)
                for key,vals in self.shufty.items():
                    if(len(vals)>1):
                        vols=[] #extract shifts
                        for val in vals:
                            vols.append(val[0])
                        if(numpy.std(vols)>1):
                            self.errors_for_combobox[pk]=1
                            diffs=[]
                            for i in range(len(vols)):
                                diff=0
                                for j in range(len(vols)):
                                    if(i!=j):
                                        diff+=vols[i]-vols[j]
                                diffs.append(diff)
                            argy=numpy.argmax(diffs)
                            spoc=vals[argy][1]
                            ii=vals[argy][2]
                            pkErr=self.peak[pk][spoc][ii]
                            self.errors_for_combobox[pk]=1



    def GetErrors(self,pk):
        ERRORS=[]
        specord='hnco','hncaco','hnca','hncoca','hncacb','hncocacb','hncocanh','hncanh','ctocsy','hcconh'
        Calpha_m1 = 0
        minus1=0
        for spec in specord:
            if(spec in self.peak[pk].keys()):
                

                if(spec=='hnco'):
                    if(len(self.peak[pk][spec])>1):
                        ERRORS.append("too many peaks in the HNCO (1 expected) %s %s" % (pk,spec))
                if(spec=='hncaco'):
                    if(len(self.peak[pk][spec])>2):
                        ERRORS.append("too many peaks in the HNCACO (2 expected) %s %s" % (pk,spec))
                if(spec=='hnca'):
                    if(len(self.peak[pk][spec])>2):
                        ERRORS.append("too many peaks in the HNCA (max 2 expected) %s %s" % (pk,spec))
                if(spec=='hncoca'):
                    if(len(self.peak[pk][spec])>1):
                        ERRORS.append("too many peaks in the HNCOCA (max 1 expected) %s %s" % (pk,spec))
                if(spec=='hncacb'):
                    if(len(self.peak[pk][spec])>4):
                        ERRORS.append("too many peaks in the HNCACB (max 4 expected) %s %s" % (pk,spec))
                    nmax=0;nmin=0;pmax=0;pmin=0;
                    posmax_shift = 0
                    for i,pk2 in enumerate(self.peak[pk][spec]):
                        if(pk2.tp=='PosMax'):
                            pmax+=1
                            posmax_shift = pk2.f3
                        if(pk2.tp=='PosMin'):
                            pmin+=1
                        if(pk2.tp=='NegMax'):
                            nmax+=1
                        if(pk2.tp=='NegMin'):
                            nmin+=1
                    if(pmin>1 or pmax>1 or nmax>1 or nmin>1):
                        ERRORS.append("MISS-ASSIGNED HNCACB: pmax %i pmin %i nmax %i nmin %i %s %s" %(pmax,pmin,nmax,nmin,pk,spec))
                    if(pmax==1 and pmin==0):
                        ERRORS.append("POSSIBLE ERROR HNCACB: no PosMin  %s %s" %(pk,spec))
                    if(nmax + nmin < 2):
                        if(posmax_shift!=0):
                            if(pk2.f3 < 47.0):
                                ERRORS.append('Note (HNCACB): Residue missing a negative CB peak and CA shift indicates glycine')

                if(spec=='hncocacb'):
                    if(len(self.peak[pk][spec])>2):
                        ERRORS.append("too many peaks in the HNCOCACB (max 2 expected) %s %s" % (pk,spec))
                    
                    for i,pk2 in enumerate(self.peak[pk][spec]):
                        if(pk2.tp=='Pos'):
                            Calpha_m1=pk2.f3
                
                if(spec=='hncocanh'):
                    if(len(self.peak[pk][spec])>2):
                        ERRORS.append("too many peaks in the HNCOCANH (max 2 expected) %s %s" % (pk,spec))

                
                if(spec=='hncanh'):
                    minus1=0
                    if(len(self.peak[pk][spec])>3):
                        ERRORS.append("too many peaks in the HNCANH (max 3 expected) %s %s" % (pk,spec))
                    
                    for i,pk_canh in enumerate(self.peak[pk][spec]):
                        if(pk_canh.tp=='minus'):
                            minus1+=1
                
                if(minus1==0):
                    if(Calpha_m1!=0 ):
                        if(Calpha_m1>60.0):
                            ERRORS.append('Potential Proline for peak i-1')


                
                

        self.GetShift(pk,tocsy=False)
        print('shufty')
        print(self.shufty)
        for key,vals in self.shufty.items(): #for each spectrum...
            if(len(vals)>1):  #for each peak...
                vols=[] #extract shifts
                for val in vals: #for each shift...
                    vols.append(val[0])
                if(numpy.std(vols)>1):
                    ERRORS.append('ERROR: inconsistent shifts: %s %.2f %.2f %s' % (key,numpy.average(vols),numpy.std(vols),vols))

                    diffs=[]
                    for i in range(len(vols)):
                        diff=0
                        for j in range(len(vols)):
                            if(i!=j):
                                diff+=vols[i]-vols[j]
                        diffs.append(diff)
                    argy=numpy.argmax(diffs)
                    spoc=vals[argy][1]
                    ii=vals[argy][2]
                    pkErr=self.peak[pk][spoc][ii]
                    ERRORS.append('%s %s %.2f %s %.2f' % (spoc,pkErr.name,pkErr.f3,pkErr.tp,pkErr.inty))

        return ERRORS

    #populate 'shufty'
    #a dictionary that holds all peak labels and checmical shifts
    def GetShift(self,pk,tocsy=True): 
        self.shufty={}
        for spec in self.peak[pk].keys():
            for i,pk3 in enumerate(self.peak[pk][spec]):
                if(tocsy==False and (spec=='ctocsy' or spec=='hcconh')):
                    continue
                lab=self.GetLab(spec,pk3.tp)
                if(lab not in self.shufty.keys()):
                    self.shufty[lab]=[]
                self.shufty[lab].append((pk3.f3p,spec,i))

    def LoadAssign(self,results):
        self.LOADASSIGN=1
        self.results=results
        self.resinv={}
        for key,val in results.items():

            #   for val in vals:
                if val not in self.resinv.keys():
                    self.resinv[val]=[]
                if(key not in self.resinv[val]):
                    self.resinv[val].append(key)
        
        print('results=%s' % self.results)
        print('resinv = %s' % self.resinv)
    

    """
    def MakeNMRstar(self):
        self.nmrStar={}
        for G1 in self.G1_nodes: #for each sequence residue in order
            resi=int(re.findall(r'\d+',G1)[0])
            resn=G1.split(str(resi))[1]
            resn3=self.p1to3[resn]

            if(G1 in self.resinv.keys()):
                if(len(self.resinv[G1])==1):
                    #print(G1,self.resinv[G1] #residue,assignment)
                    pk=self.resinv[G1][0]    #nmrID
                    if(resi-1 in self.seq.keys()): #is there a preceeding residue?
                        if(self.seq[resi-1]=='P' or resi-1==1): #if the previous residue was a proline, or start of chain write i-1s
                            #print(G1,resi,resn,self.seq[resi-1],self.seq[resi])
                            resnP=self.seq[resi-1]
                            resn3P=self.p1to3[resnP]
                            for spec in self.peak[pk].keys():
                                for i,pk3 in enumerate(self.peak[pk][spec]):
                                    if(resi-1 not in self.nmrStar.keys()):
                                        self.nmrStar[resi-1]={}
                                    lab=self.GetLab(spec,pk3.tp)
                                    if(lab=='CO(i-1)' and spec=='hnco'):
                                        self.nmrStar[resi-1]['CO']=pk3.f3
                                    elif(lab=='CA(i-1)' and spec=='hnca'):
                                        self.nmrStar[resi-1]['CA']=pk3.f3
                                    elif(lab=='CB(i-1)' and spec=='hncacb'):
                                        self.nmrStar[resi-1]['CB']=pk3.f3

                    if(resi not in self.nmrStar.keys()):
                        self.nmrStar[resi]={}

                    if(self.refSpec in self.peak[pk].keys()):
                        self.nmrStar[resi]['N']=self.peak[pk][self.refSpec][0].f2
                        self.nmrStar[resi]['HN']=self.peak[pk][self.refSpec][0].f1
                    for spec in self.peak[pk].keys():
                        for i,pk3 in enumerate(self.peak[pk][spec]):
                            lab=self.GetLab(spec,pk3.tp)
                            if(lab=='CO(i)'):
                                self.nmrStar[resi]['CO']=pk3.f3p
                            elif(lab=='CA(i)' and spec=='hnca'):
                                self.nmrStar[resi]['CA']=pk3.f3p
                            elif(lab=='CB(i)'):
                                self.nmrStar[resi]['CB']=pk3.f3p
    """

    #write peak list files as output
    def assemble(self,outfile,results):
        self.LoadAssign(results) #make resultsinv dictionary
        self.MakeNMRstar()  #assemble nmrStar dictionary

        self.RenamePeaks()  #rename peak array


        if(os.path.exists('peaklists')==0):
            os.system('mkdir peaklists')
        outfile=os.path.join('peaklists',outfile)
        print(outfile)
        cnt=1 #print(BMRB type summary)
        outy=open(outfile,'w')
        for resi in self.nmrStar.keys():
            resn=self.seq[resi]
            resn3=self.p1to3[resn]
            for atom in ('HN','N','CO','CA','CB'):
                try:
                    outy.write('%i %i %s %s %s %f\n' % (cnt,resi,resn,resn3,atom,self.nmrStar[resi][atom]));cnt+=1
                except:
                    pass
        outy.close()


        #print(talos tab type summary)
        outy=open(outfile+'.tab','w')
        # outy.write('DATA FIRST_RESID %i\n' %  self.FirstResidue)
        outy.write('DATA SEQUENCE ')
        #for resi in self.resi:
        #    outy.write('%s' % (self.seq[resi]))

        seq=''
        mapp={}
        cnt=1
        for resi in self.nmrStar.keys():
            if(resi in self.nmrStar.keys()):
                mapp[resi]=cnt
                seq+=self.seq[resi]
                cnt+=1
        outy.write('%s' % (seq))



        outy.write('\n')
        outy.write('VARS   RESID RESNAME ATOMNAME SHIFT\n')
        outy.write('FORMAT %s4d   %s1s     %s4s      %s8.3f\n' % ('%','%','%','%'))
        #form='%s4d   %s1s     %s4s      %s8.3f\n'
       
        for resi in self.nmrStar.keys():
            resn=self.seq[resi]
            resn3=self.p1to3[resn]
            if(resi not in self.nmrStar.keys()):
                continue
            for atom in self.nmrStar[resi]:
                outy.write('%4d   %1s     %4s      %8.3f\n' % (mapp[resi],self.seq[resi],atom,self.nmrStar[resi][atom]))
        outy.close()

        #print(HN)
        outy=open(outfile+'.HN','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    outy.write('%i%s\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['N'],self.nmrStar[resi]['HN']))
        outy.close()


        def FindTyp(typ,vals):
            for koi,vol in vals.items():
                if(len(typ)==1): #if type has just one letter..
                    if(koi==typ):
                        return vol
                elif(len(typ)>=2):
                    if(koi==typ):
                        return vol
            return False
        #print(HC)
        outy=open(outfile+'.CH','w')
        for resi,vals in self.nmrStar.items():
            for t in 'A','B','G','D','E':
                for b in '','1','2','3':
                    C=FindTyp('C'+t+b,vals)  #get the C chemical shift associated with Ctb
                    if(C!=False):
                        for a in '','1','2','3':
                            H=FindTyp('H'+t+b+a,vals)  #seach the lists to find H+t+b+a, the matching proton.
                            if(H!=False):
                                outy.write('%s\t%s\t%s\n' % (str(resi)+self.seq[resi]+'C'+t+b+'-'+'H'+t+b+a,C,H))


            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    outy.write('%i%sH-C\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N']))
        outy.close()

        #print(hnco)
        outy=open(outfile+'.hnco','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    if(resi-1 in self.nmrStar.keys()):
                        if('CO' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CO\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CO']))
        outy.close()

        #print(hncaco)
        outy=open(outfile+'.hncaco','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    if(resi-1 in self.nmrStar.keys()):
                        if('CO' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CO\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CO']))
                    if('CO' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CO\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CO']))
        outy.close()

        #print(hnca)
        outy=open(outfile+'.hnca','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    if('CA' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CA']))
                    if(resi-1 in self.nmrStar.keys()):
                        if('CA' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
        outy.close()

        #print(hncoca)
        outy=open(outfile+'.hncoca','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    #if('CA' in self.nmrStar[resi].keys()):
                    #    outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CA']))
                    if(resi-1 in self.nmrStar.keys()):
                        if('CA' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
        outy.close()

        #print(hncacb)
        outy=open(outfile+'.hncacb','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    if('CA' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CA']))
                    if('CB' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CB\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CB']))

                    if(resi-1 in self.nmrStar.keys()):
                        if('CA' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
                        if('CB' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CB\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
        outy.close()

        #print(hncocacb)
        outy=open(outfile+'.hncocacb','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    if('CA' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CA']))
                    if('CB' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CB\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CB']))
                    if(resi-1 in self.nmrStar.keys()):
                        if('CA' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
                        if('CB' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CB\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
        outy.close()






    def LoadAssignment(self,results):
        self.LOADASSIGN=1
        self.results=results
        self.resinv={}
        
        for key,val in results.items():
            
            if val not in self.resinv.keys():
                self.resinv[val]=[]
            if(key not in self.resinv[val]):
                self.resinv[val].append(key)
        
        print('results=%s' % self.results)
        print('resinv = %s' % self.resinv)
    
    def StarToBMRB(self,typ):
            tt=typ  #convert nmrStar to BMRB
            if(tt=='CO'):
                tt='C'
            if(tt=='HN'):
                tt='H'
            if(typ=='HBD'):
                return 'HB'
            if(typ=='CBG'):
                return 'CB'
            
            return tt


    def TestAlias(self,spec,resi,typ,ppm,dim):
        
        tt=self.StarToBMRB(typ)

        tig=0
        if(resi in self.shiftx2.keys()):
            for a in '','1','2','3':  #try some combinations if the first option is not there.
                if(tt+a in self.shiftx2[resi][tt[0]].keys()):
                    shift=self.shiftx2[resi][tt[0]][tt+a]
                    tig=1
                    break

        else:
            
            #print(typ,tt,len(typ))

            #print(self.bmrb[typ[0]][tt])
            for resn,shift,std in self.bmrb[typ[0]][tt]:
                #open aprint(resn,shift,std)
                if(self.p1to3[self.seq[resi]]==resn):
                    tig=1
                    #print (shift)
                    break
        #sys.exit(100)
        if(tig==0):
            print('ALIAS FAIL:',spec,resi,typ,ppm,dim)
            return ppm
        
        
        dd,dmax,dmin,vals=self.spec[spec].GetAliasRange(dim) #for specified dimension, get alias values
        
        #get some permuatations of aliasing...
        ppmUp=ppm+(dmax-dmin+dd)
        ppmDown=ppm-(dmax-dmin+dd)
        ppmUp2=ppmUp+(dmax-dmin+dd)
        ppmDown2=ppmDown-(dmax-dmin+dd)

        ppmTest=ppm,ppmUp,ppmDown,ppmUp2,ppmDown2
        ppmTest=numpy.array(ppmTest)
        test=numpy.fabs(ppmTest-shift) #sutract expectation value
        argy=numpy.argmin(test)  #find the closest....
        if(argy!=0): #if the closest needs aliasing...
            print('NMRSTAR ALIASING:',spec,resi,typ,self.seq[resi],ppm,ppmTest[argy])
            #print(ppmTest)
            #print(test)
            #print(argy)
        #if(resi==86 and typ=='N'):
        #    print('NMRSTAR ALIASING:',spec,resi,typ,self.seq[resi],ppm,ppmTest[argy])
        #    #sys.exit(100)
        return ppmTest[argy]  #return adjusted shift.


        


    def MakeNMRstar(self):
        print('Creating NMRstar dict...')
        self.nmrStar={}
        #print(self.G1_nodes,len(self.G1_nodes))
        for G1 in self.G1_nodes: #for each sequence residue (atom) in order
            resi=int(G1[:-1])
            #resi=int(re.findall(r'\d+',G1)[0])
            resn=G1[-1]  #residue single letter
            resn3=self.p1to3[resn]  #residue 3 letter
            #print(G1,self.resinv.keys())
            if(G1 in self.resinv.keys()):
                if(len(self.resinv[G1])==1):
                    #print(G1,self.resinv[G1] #residue,assignment)
                    pk=self.resinv[G1][0]    #nmrID
                    if( (resi-1) in self.seq.keys()): #is there a preceeding residue?
                        #print(G1,pk,resi-1)
                        #if(self.seq[resi-1]=='P' or resi-1==1): #if the previous residue was a proline, or start of chain write i-1s
                        #print(G1,resi,resn,self.seq[resi-1],self.seq[resi])
                        resnP=self.seq[resi-1]
                        resn3P=self.p1to3[resnP]
                        for spec in self.peak[pk].keys(): #for all associated spectra...    
                            for i,pk3 in enumerate(self.peak[pk][spec]): #go over all peaks
                                lab=self.GetLab(spec,pk3.tp)   #get the expected label.
                                #print(G1,pk,lab,spec,pk3.name,pk3.f3p)
                                if(lab==''):  #peak is not labelled.
                                    continue
                                if('(i-1)' not in lab):
                                    continue
                                if( (resi-1) not in self.nmrStar.keys()):
                                    self.nmrStar[resi-1]={}
                
                                typ=lab.split('(i-1)')[0]

                                if(spec=='ctocsy'):  #skip this one as it is a repeat of HNCA.
                                    if('CA' in lab):
                                        continue
                            
                                self.nmrStar[resi-1][typ]=self.TestAlias(spec,resi-1,typ,pk3.f3,0)

                    if(resi not in self.nmrStar.keys()):
                        self.nmrStar[resi]={}

                    if(pk not in self.peak.keys()):
                        continue
                    if(self.refSpec in self.peak[pk].keys()):
                        #do the HN
                        self.nmrStar[resi]['N']=self.TestAlias(self.refSpec,resi,'N',self.peak[pk][self.refSpec][0].f2,1)
                        self.nmrStar[resi]['HN']=self.TestAlias(self.refSpec,resi,'HN',self.peak[pk][self.refSpec][0].f1,2)

                    for spec in self.peak[pk].keys():
                        for i,pk3 in enumerate(self.peak[pk][spec]):
                            lab=self.GetLab(spec,pk3.tp)

                            if(lab==''):  #peak is not labelled.
                                continue
                            if(lab[0]!='C'):  #only know residue i carbons
                                continue
                            if('(i)' not in lab):
                                continue
                            if( (resi) not in self.nmrStar.keys()):
                                self.nmrStar[resi]={}
                            typ=lab.split('(i)')[0]

                            self.nmrStar[resi][typ]=self.TestAlias(spec,resi,typ,pk3.f3,0)

                            #if(lab=='CO(i)'):
                            #    self.nmrStar[resi]['CO']=self.TestAlias(spec,resi,'CO',pk3.f3,0)
                            #elif(lab=='CA(i)' and spec=='hnca'):
                            #    self.nmrStar[resi]['CA']=self.TestAlias(spec,resi,'CA',pk3.f3,0)
                            #elif(lab=='CB(i)'):
                            #    self.nmrStar[resi]['CB']=self.TestAlias(spec,resi,'CB',pk3.f3,0)
        



        #for key,vals in self.nmrStar.items():
        #    for koi,vols in vals.items():
        #        print(key,koi,vols)
        #print(self.ShowPeak(self.peak['75']['ctocsy']))
        #sys.exit(100)
        #print('bb',self.nmrStar[85]['N'])
        #sys.exit(100)

    #write peak list files as output
    def assemble_assignment(self,outfile,results):
        self.LoadAssignment(results) #make resultsinv dictionary
        self.MakeNMRstar()  #assemble nmrStar dictionary

        self.RenamePeaks()  #rename peak array


        if(os.path.exists('peaklists')==0):
            os.system('mkdir peaklists')
        outfile=os.path.join('peaklists',outfile)
        print(outfile)
        cnt=1 #print(BMRB type summary)
        outy=open(outfile,'w')
        for resi in self.nmrStar.keys():
            resn=self.seq[resi]
            resn3=self.p1to3[resn]
            for atom in ('HN','N','CO','CA','CB'):
                try:
                    outy.write('%i %i %s %s %s %f\n' % (cnt,resi,resn,resn3,atom,self.nmrStar[resi][atom]));cnt+=1
                except:
                    pass
        outy.close()


        #print(talos tab type summary)
        outy=open(outfile+'.tab','w')
        outy_ca = open(outfile+'.ca','w')
        outy_cb = open(outfile+'.cb','w')

        # outy.write('DATA FIRST_RESID %i\n' %  self.FirstResidue)
        outy.write('DATA SEQUENCE ')
        for resi in self.resi:
            outy.write('%s' % (self.seq[resi]))
        outy.write('\n')
        outy.write('VARS   RESID RESNAME ATOMNAME SHIFT\n')
        outy.write('FORMAT %s4d   %s1s     %s4s      %s8.3f\n' % ('%','%','%','%'))
        form='%s4d   %s1s     %s4s      %s8.3f\n'
        for resi in self.nmrStar.keys():
            resn=self.seq[resi]
            resn3=self.p1to3[resn]
            for atom in ('HN','N','CO','CA','CB'):

                try:
                    outy.write('%4d   %1s     %4s      %8.3f\n' % (resi,self.seq[resi],atom,self.nmrStar[resi][atom]))
                except:
                    pass
        outy.close()
        
        for resi in self.nmrStar.keys():
            resn=self.seq[resi]
            resn3=self.p1to3[resn]

            try:
                outy_ca.write(resi+' '+self.nmrStar[resi]['CA'])
            except:
                pass

            try:
                outy_cb.write(resi+' '+self.nmrStar[resi][atom])
            except:
                pass


        
        
        
        
        
        
        outy_ca.close()
        outy_cb.close()




        #print(HN)
        outy=open(outfile+'.HN','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    outy.write('%i%sHN-N\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N']))
        outy.close()

        #print(hnco)
        outy=open(outfile+'.hnco','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    if(resi-1 in self.nmrStar.keys()):
                        if('CO' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CO\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CO']))
        outy.close()

        #print(hncaco)
        outy=open(outfile+'.hncaco','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    if(resi-1 in self.nmrStar.keys()):
                        if('CO' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CO\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CO']))
                    if('CO' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CO\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CO']))
        outy.close()

        #print(hnca)
        outy=open(outfile+'.hnca','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    if('CA' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CA']))
                    if(resi-1 in self.nmrStar.keys()):
                        if('CA' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
        outy.close()

        #print(hncoca)
        outy=open(outfile+'.hncoca','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    #if('CA' in self.nmrStar[resi].keys()):
                    #    outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CA']))
                    if(resi-1 in self.nmrStar.keys()):
                        if('CA' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
        outy.close()

        #print(hncacb)
        outy=open(outfile+'.hncacb','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    if('CA' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CA']))
                    if('CB' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CB\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CB']))

                    if(resi-1 in self.nmrStar.keys()):
                        if('CA' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
                        if('CB' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CB\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
        outy.close()

        #print(hncocacb)
        outy=open(outfile+'.hncocacb','w')
        for resi in self.nmrStar.keys():
            if('HN' in self.nmrStar[resi].keys()):
                if('N' in self.nmrStar[resi].keys()):
                    if('CA' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CA']))
                    if('CB' in self.nmrStar[resi].keys()):
                        outy.write('%i%sHN-N-CB\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi]['CB']))
                    if(resi-1 in self.nmrStar.keys()):
                        if('CA' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CA\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
                        if('CB' in self.nmrStar[resi-1].keys()):
                            outy.write('%i%sHN-N-CB\t%f\t%f\t%f\n' % (resi,self.seq[resi],self.nmrStar[resi]['HN'],self.nmrStar[resi]['N'],self.nmrStar[resi-1]['CA']))
        outy.close()















    def AnalPeak(self,pk):
        self.ShiftScore_i(pk)
        resns,probs=self.ResidueProb()
        self.ShiftScore_im1(pk)
        resns2,probs2=self.ResidueProb(resns=resns,p=True)
        return resns,probs,probs2


    def ShiftScore_i(self,pk):
        self.scr={}
        #print('keys = %s' %self.peak[pk].keys())
        if(self.refSpec in self.peak[pk].keys()):
            #print('refspec')
            for i,pk3 in enumerate(self.peak[pk][self.refSpec]):
                # print(pk3.name,pk3.f2,pk3.f1)
                self.DoScr('N','N',pk3.f2)
                self.DoScr('H','H',pk3.f1+self.Hfudge)  #NOTE! REFERENCING FUDGE!
                #self.DoScr('H','H',pk3.f1)  #NOTE! REFERENCING FUDGE!
                break
        else: #otherwise default to hnco
            key = list(self.peak[pk].keys())[0]
            for i,pk3 in enumerate(self.peak[pk][key]):
                # print(pk3.name,pk3.f2,pk3.f1)
                self.DoScr('N','N',pk3.f2)
                #self.DoScr('H','H',pk3.f1+1)  #NOTE! REFERENCING FUDGE!
                self.DoScr('H','H',pk3.f1+self.Hfudge)  #NOTE! REFERENCING FUDGE!
                break

        if('hncaco' in self.peak[pk].keys()):
            for i,pk3 in enumerate(self.peak[pk]['hncaco']):
                if(pk3.tp!='main'):
                    self.DoScr('C','C',pk3.f3p)
                    break

        if('hnca' in self.peak[pk].keys()):
            for i,pk3 in enumerate(self.peak[pk]['hnca']):
                if pk3.tp=='main':
                    self.DoScr('C','CA',pk3.f3p)
                    break

        if('hncacb' in self.peak[pk].keys()):
            for i,pk3 in enumerate(self.peak[pk]['hncacb']):
                if pk3.tp=='NegMax':
                    self.DoScr('C','CB',pk3.f3p)
                    break

        # if('hncocacb' in self.peak[pk].keys()):
        #     for i,pk3 in enumerate(self.peak[pk]['hncocacb']):
        #         if pk3.tp=='Neg':
        #             self.DoScr('C','CB',pk3.f3p)
                    # break
        # if('hncanh' in self.peak[pk].keys()):
        #     for i,pk3 in enumerate(self.peak[pk]['hncanh']):
        #         if pk3.tp=='diag':
        #             self.DoScr('N','N',pk3.f2)
        #             self.DoScr('H', 'H', pk3.f3p)
        #             break
        # if('hncocanh' in self.peak[pk].keys()):
        #     for i,pk3 in enumerate(self.peak[pk]['hncocanh']):
        #         if pk3.tp=='diag':
        #             self.DoScr('N','N',pk3.f2)
        #             self.DoScr('H', 'H', pk3.f3p)
        #             break


    def ShiftScore_im1(self,pk):
        self.scr={}
        if('hnco' in self.peak[pk].keys()):
            for i,pk3 in enumerate(self.peak[pk]['hnco']):
                self.DoScr('C','C',pk3.f3p)
                break

        if('hnca' in self.peak[pk].keys()):
            for i,pk3 in enumerate(self.peak[pk]['hnca']):
                if pk3.tp=='': #CA i-1
                    self.DoScr('C','CA',pk3.f3p)
                    break

        if('hncacb' in self.peak[pk].keys()):
            for i,pk3 in enumerate(self.peak[pk]['hncacb']):
                if pk3.tp=='NegMin': #CB i-1
                    self.DoScr('C','CB',pk3.f3p)
                    break

        if('cbcaconh' in self.peak[pk].keys()):
            for i,pk3 in enumerate(self.peak[pk]['cbcaconh']):
                # print(pk3)
                if pk3.tp!='high': #CB i-1
                    self.DoScr('C','CB',pk3.f3p)
                    break

        #print('scr i-1:')
        #for key,val in self.scr.items():
        #    print(key,val)
        #    print(self.bmrbC[self.p1to3[key]])      

        if('ctocsy' in self.peak[pk].keys()): #loop over all residues and score best matching assignment
            self.DoScrTOCSY('C',pk,self.bmrbC)
            
        if('hcconh' in self.peak[pk].keys()):#loop over all residues and score best matching assignment
            self.DoScrTOCSY('H',pk,self.bmrbH)
            #print('done this')

        # print('dfdfafdafda')
        #print('scr i-1:')
        #for key,val in self.scr.items():
        #    print(key,val)
        #    try:
        #        print(self.bmrbC[self.p1to3[key]])            
        #    except:
        #        pass

    def PeakShifts(self):
        print('-------------------------------------------------')
        print('Analysing assignment options from chemical shifts')
        # try:
        #     print(self.peak['3H-N'])
        # except:
        #     print(self.peak.keys())
        #     exit()
        

        #self.SetBMRB()
        keys=self.peak.keys()
        keys=sorted(keys)



        self.shift={}
        self.shift2={}


        for pk in keys:
            print('---------------------------------------')
            print('Analysing peak:',pk)
            self.ShiftScore_i(pk)
            #print('self.scr.keys = %s' % self.scr.keys())
            resns,probs=self.ResidueProb()
            if len(probs) == 0:
                print('No assignment options for',pk)
            # print('shift dict before')
            # print(self.shift)
            self.ResidueClassify(pk,self.shift,resns,probs)
            
            #print(self.shift[pk])
            #print(self.shift[pk][:,1])
            #argy=numpy.argsort(self.shift[pk][:,1])
            #print(self.shift[pk])
            print('i:',self.shift[pk][:,0])
            # print('shift dict before')
            # print(self.shift)
            
            #
            # if pk == "43H-N":
            #     print(self.shift["43H-N"])
            #     print(self.tolMin)
            #     #print(self.shift2["43H-N"])
            #     exit()

            self.ShiftScore_im1(pk)

            #print('scr i-1:')
            #for key,val in self.scr.items():
            #    print(key,val)
            #    print(self.bmrbC[self.p1to3[key]])            


            resns,probs=self.ResidueProb(resns=resns,p=True)
            self.ResidueClassify(pk,self.shift2,resns,probs,p=True)
            #print(self.shift2[pk])
            print('i-1:',self.shift2[pk][:,0])
        # try:
        #     print('self.shift dictionary')
        #     print(self.shift['3H-N'])
        # except:
        #     print('This is the problem!!')
        #     sys.exit(100)

            #cheat functions: make sure known residue type has not been excluded
            #self.CheckErrorShift(pk,self.shift)
            #self.CheckErrorShift(pk,self.shift2,im1=True)
        


    def ResidueProb(self,resns=[],p=False):
        #print('Calculating residue probabilities based on BMRB')
        #go through possible residues and work out relative probabilities
        #based on BMRB statistics.
        if(len(resns)==0):
   
            resns=self.scr.keys() #get rid of prolines
            #print('resns')
            #print(resns)
            resns = list(resns)
            sorted(resns)


        probs=[]
        for resn in resns:
            #print(resn)
            prob=1. #start at 1...
            if(p==False): #set p prob to zero if needed
                if(resn=='P'):
                    prob=0
            if(resn in self.scr.keys()):
                for tp in self.scr[resn]: #for all entries in dictionary
#                    print(resn,tp,self.scr[resn][tp])
                    # if resn == "43H-N":
                       # print(tp)
       # print(self.scr[resn][tp])

                    prob*=self.scr[resn][tp]
            #if resn == "7H-N":
            #exit()
            probs.append(prob)
        # print(resns[10], probs[10])
        return resns,probs


    def ResidueClassify(self,pk,shift,resns,probs,p=False):
        #argprob=numpy.argmax(probs)
        #probmax=probs[argprob]

        
        #print(probs, pk)
        shift[pk] = []
        
        argprob=numpy.argsort(probs)  #sort probabilities (return arg)
        #print(argprob)
        probs=numpy.array(probs)
        #print(probs[argprob])
        #print(argprob)
        probmax=probs[argprob[-1]]    #get max prob for all residues
        #print(probs[argprob[-1]]/probs[argprob[-2]])

        for i in range(3): #RULE 1: if there is a sudden drop off, we are certian of identification. Go to depth=3
            self.tolSin=20
            if(probs[argprob[-1-i]]/probs[argprob[-2-i]]>self.tolSin):
                resn=resns[argprob[-1-i]]
                for j in range(i+1):
                    reson=resns[argprob[-1-j]]
                    if(p):
                        shift[pk].append((reson,probs[argprob[-1-j]]))
                    else:
                        if(reson!='P'):
                            shift[pk].append((reson,probs[argprob[-1-j]]))
                shift[pk]=numpy.array(shift[pk])
                return

        if(probmax>self.tolMax):#RULE 2: if high probs are present, take a bunch above tolmin
            for ii in range(len(argprob)):
                i=argprob[-ii-1]
                resn=resns[i]
                #for i,resn in enumerate(resns):
                if(probs[i]>self.tolMin):
                    #shift[pk].append(self.p3to1[resn])
                    if(p):
                        shift[pk].append((resn,probs[i]))
                    else:
                        if(resn!='P'):
                            shift[pk].append((resn,probs[i]))
            shift[pk]=numpy.array(shift[pk])
            return

        for ii in range(len(argprob)):
            i=argprob[-ii-1]
            resn=resns[i]
            #print(resn,probs[i])
            #for i,resn in enumerate(resns): #RULE 3: we've no clue. take everybody
            #shift[pk].append(self.p3to1[resn])
            if(p):
                shift[pk].append((resn,probs[i]))
            else:
                if(resn!='P'):
                    shift[pk].append((resn,probs[i]))
        #print (shift[pk])
        shift[pk]=numpy.array(shift[pk])


    """
    #Cheat function - make sure residue type has not been excluded by shift filter
    def CheckErrorShift(self,pk,shift,im1=False):
        resi=int(re.findall(r'\d+',pk)[0])
        print(resi)
        if(resi==21):
            resi=28
        elif(resi==28):
            resi=21
        elif(resi==72):
            resi=31
        elif(resi==31):
            resi=72
        elif(resi==73):
            resi=69
        elif(resi==69):
            resi=73

        if(im1): #subtract 1 to get i-1
            resi-=1

        cheat=self.seq[resi]
        #print(pk,resi,cheat)
        if(cheat not in shift[pk]):
            print('shit - fucked up chemical shifts')
            if(im1):
                print('im1')
            print(shift[pk])
            print(cheat)
            #if(pk!='53H-N'):
            #    sys.exit(100)
    """

    def SetBMRB(self):

        
        import pathlib
        #print()
        bmrbfile=os.path.join(pathlib.Path(__file__).parent.resolve(),'bmrb.txt')
        #print(bmrbfile)
        if(os.path.exists(bmrbfile)==False):
            print('Cannot find BMRB file. Check installation.')
            sys.exit(100)
        #sys.exit(100)

        self.bmrb={}
        self.bmrbC={}
        self.bmrbH={}
        #print(os.getcwd())
        #print(os.path.realpath(os.path.dirname(__file__)))

        inny=open(bmrbfile)
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                resn=test[0]
                typ=test[1]
                atom=test[2]
                shift=float(test[6])
                std=float(test[7])

                if(atom not in self.bmrb.keys()):
                    self.bmrb[atom]={}
                if(typ not in self.bmrb[atom]):
                    self.bmrb[atom][typ]=[]
                self.bmrb[atom][typ].append((resn,shift,std))

                if(atom=='C'):  ## CB enabling COs in the bmrbC for plotting in the spinsystem tab
                # if(atom=='C' and shift < 100):
                    if(resn not in self.bmrbC.keys()):
                        self.bmrbC[resn]={}

                    self.bmrbC[resn][typ]=shift,std
                if(atom=='H'):
                    if(resn not in self.bmrbH.keys()):
                        self.bmrbH[resn]={}
                    self.bmrbH[resn][typ]=shift,std
                    #if(resn=='VAL' and typ=='CG'):
                    #    self.bmrbC[resn]['Cg']=shift,std
                    #if(resn=='LEU' and typ=='CD'):
                    #    self.bmrbC[resn]['Cd']=shift,std

    def GetProbShift(self,shiftxresi,resD,probsA,probs2A):
                resiI=self.shiftx2[shiftxresi]['resn']
                resiIm1=self.shiftx2[shiftxresi-1]['resn']

                if(resiI not in resD.keys()):
                    return 0
                if(resiIm1 not in resD.keys()):
                    return 0
      
                probI=probsA[resD[resiI]]
                probIm1=probs2A[resD[resiIm1]]
                prob=probI*probIm1
                return prob

    def GetProbShiftX2(self,pk,shiftxresi,cerr):
            prob=1.
            for lab in self.shufty.keys(): #CO(i),CO(i-1)...etc
                ref=self.shufty[lab][0][0]  #chemical shift of residue associated with this.
                if(lab=='CO(i)'):
                    if('C' in self.shiftx2[shiftxresi].keys()):
                        if('CO' in self.shiftx2[shiftxresi]['C'].keys()):
                            ppm=self.shiftx2[shiftxresi]['C']['CO']
                            prob*=numpy.exp(-(ppm-ref)**2./(2.*cerr**2.))
                elif(lab=='CA(i)'):
                    if('C' in self.shiftx2[shiftxresi].keys()):
                        if('CA' in self.shiftx2[shiftxresi]['C'].keys()):
                            ppm=self.shiftx2[shiftxresi]['C']['CA']
                            prob*=numpy.exp(-(ppm-ref)**2./(2.*cerr**2.))

                elif(lab=='CB(i)'):
                    if('C' in self.shiftx2[shiftxresi].keys()):
                        if('CB' in self.shiftx2[shiftxresi]['C'].keys()):
                            ppm=self.shiftx2[shiftxresi]['C']['CB']
                            prob*=numpy.exp(-(ppm-ref)**2./(2.*cerr**2.))

                elif(lab=='CO(i-1)'):
                    if(shiftxresi-1 in self.shiftxresis):
                        if('C' in self.shiftx2[shiftxresi-1].keys()):
                            if('CO' in self.shiftx2[shiftxresi-1]['C'].keys()):
                                ppm=self.shiftx2[shiftxresi-1]['C']['CO']
                                prob*=numpy.exp(-(ppm-ref)**2./(2.*cerr**2.))
                elif(lab=='CA(i-1)'):
                    if(shiftxresi-1 in self.shiftxresis):
                        if('C' in self.shiftx2[shiftxresi-1].keys()):
                            if('CA' in self.shiftx2[shiftxresi-1]['C'].keys()):
                                ppm=self.shiftx2[shiftxresi-1]['C']['CA']
                                prob*=numpy.exp(-(ppm-ref)**2./(2.*cerr**2.))
                

                elif(lab=='CB(i-1)'):
                    if(shiftxresi-1 in self.shiftxresis):
                        if('C' in self.shiftx2[shiftxresi-1].keys()):
                            if('CB' in self.shiftx2[shiftxresi-1]['C'].keys()):
                                ppm=self.shiftx2[shiftxresi-1]['C']['CB']
                                prob*=numpy.exp(-(ppm-ref)**2./(2.*cerr**2.))
                              
                """
                elif('ctocsy' in self.spec.keys()):
                    if(shiftxresi-1 in self.shiftxresis):
                        if('C' in self.shiftx2[shiftxresi-1].keys()):
                            if('CB' in self.shiftx2[shiftxresi-1]['C'].keys()):
                                
                                if('ctocsy' not in self.peak[pk].keys()):
                                    continue
                                tig=0
                                for pk3 in self.peak[pk]['ctocsy']:
                                    if(pk3.tp=='CA(i-1)'):
                                        tig=1
                                        break
                                if(tig==0):
                                    continue

                                #residue type
                                resnT=self.shiftx2[shiftxresi-1]['resn']
                                #map TOCSY peaks onto this residue and classify.
                                resT=self.DoScrTOCSY('C',pk,self.shiftx2[shiftxresi-1]['C'],cerr=cerr,resn=resnT)
                                #resT=self.MatchTOCSY('C',self.p1to3[resnT],pk)
                                print (shiftxresi,pk,resT)
                                if('CB' not in resT.keys()):
                                    #prob*=0
                                    continue
                                else:
                                    refT=resT['CB'][0]
                                    ppmT=self.shiftx2[shiftxresi-1]['C']['CB']
                                    prob*=numpy.exp(-(ppmT-refT)**2./(2.*cerr**2.))
                
                                #print(self.shiftx2[shiftxresi-1])
                                #sys.exit(100)
                                #refT=XXXXX
                                #ppm=self.shiftx2[shiftxresi-1]['C']['CB']
                                #prob*=numpy.exp(-(ppm-refT)**2./(2.*cerr**2.))
                                #take the one that comes out as CB, and calculate the CB probability.
                """
            return prob

    #compare a peak list to either a BMRB file from a template
    #or to the BMRB overall statistics.
    def CompareShiftx2(self,pk):


        self.ReadShiftx2()  #read either template or BMRB file.


        #not used.
        self.GetShift(pk) #assembly experimental shufty array (no TOCSY)

        #go through peak lists directly for i/i-1 prob dists o
        resnsA,probsA,probs2A=self.AnalPeak(pk) 
        resD={}
        for i,res in enumerate(resnsA):
            resD[res]=i


        probs=[]
        resnsAll=[]
        prob1=[]
        prob2=[]

        if(len(self.template)==0): #increae widths if using BMRB only.
            cerr=3
        else: #otherwise use this.
            cerr=2

        self.USETOCSY=False
        #print("HELLO")
        #print(self.shufty)
        #print(self.shiftx2)
        print (self.shiftxresis)
        for shiftxresi in self.shiftxresis: #for each residue in the database
            if(shiftxresi==self.FirstResidue): #we cannot detect the 1st residue
                continue
            resnsAll.append(str(shiftxresi)+self.shiftx2[shiftxresi]['resn']) #append residue label.

            #print(shiftxresi,self.FirstResidue)
            prob1.append(self.GetProbShiftX2(pk,shiftxresi,cerr))#use reference files
            if(self.USETOCSY): #combine generalised residue probabilities
                try: 
                    prob2.append(self.GetProbShift(shiftxresi,resD,probsA,probs2A)) #include TOCSY
                except:
                    prob2.append(0)
            else:
                prob2.append(0)

        prob1=numpy.array(prob1) #shiftx2
        prob2=numpy.array(prob2) #full list
        resnsAll=numpy.array(resnsAll)
        #1. try using just the tocsy.
        thresh=0.05
        while(1==1):
            print(thresh)
            probs=prob2[prob2>thresh]
            resns=resnsAll[prob2>thresh]
            if(len(probs)==0):
                probs=prob1[prob1>thresh]
                resns=resnsAll[prob1>thresh]
            if(len(probs)!=0):
                break
            thresh/=5
 
        print('resn',resns)
        print('prob',probs)

        #resns=numpy.array(resns)
        #probs=numpy.array(probs)
        argy=numpy.flip(numpy.argsort(probs))
        probs=probs[argy]
        resns=resns[argy]
        
        return resns,probs


    def ReadShiftx2(self):
        print('Getting reference shifts',self.template)
        self.shiftx2={}
        if(self.template!='' and os.path.exists(self.template)):
            inny=open(self.template)
            for line in inny.readlines():
                test=line.split()
                if(len(test)>0):
                    resi=int(test[1])
                    resn=test[2]
                    atom=test[3][0]
                    typ=test[3]
                    ppm=float(test[4])
                    if(resn!='PRO'):
                        if(resi not in self.shiftx2.keys()):
                            self.shiftx2[resi]={}
                            self.shiftx2[resi]['resn']=self.p3to1[resn]
                        if(atom not in self.shiftx2[resi].keys()):
                            self.shiftx2[resi][atom]={}
                        self.shiftx2[resi][atom][typ]=ppm
            inny.close()

            ######For Gly, we are getting 3 entries for proton for some reason, 1/2 and average.
            #remove average.
            for resi,atoms in self.shiftx2.items():
                if(self.shiftx2[resi]['resn']!='G'):
                    continue
                for a,typs in atoms.items():
                    if(a!='H'):
                        continue
                    cnt=0
                    for typ,shift in typs.items():
                        if(len(typ)>=2):
                            if(typ[1]=='A'):
                                cnt+=1
                    if(cnt==3):
                        del self.shiftx2[resi]['H']['HA']

        else: #load in the BMRB as a reference
            print('Using BMRB as reference')
            #not using TOCSY?
            #try:
            #    self.resi
            #except:
            #    self.shiftxresis=[]
            #    return

            for resi in self.resi:
                resn=self.seq[resi]
                for atom in 'C','H':
                    if(atom in self.bmrb.keys()):
                        for typ in ('CA','C','CB','H','N'):
                            if typ in self.bmrb[atom].keys():
                                if(resi not in self.shiftx2.keys()):
                                    self.shiftx2[resi]={}
                                    self.shiftx2[resi]['resn']=resn
                                if(atom not in self.shiftx2[resi].keys()):
                                    self.shiftx2[resi][atom]={}

                                for (resnb,shiftb,stdb) in self.bmrb[atom][typ]:
                                    if(self.p3to1[resnb]==resn):
                                        #print(resi,resn,shiftb,self.p3to1[resn])
                                        self.shiftx2[resi][atom][typ]=shiftb

        self.shiftxresis=self.shiftx2.keys()
        self.shiftxresis=sorted(self.shiftxresis)
