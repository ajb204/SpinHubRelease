#!/opt/local/bin/python

###################################################################
# Deconvolve nmr spectrum
###################################################################

import numpy,sys,os,string,math,wx
import nmrglue as ng
#from baldwinStd import readfile



def ParseFlt(infile,param):
    inny=open(infile,'r')
    for line in inny.readlines():
        test=line.split()
        if(len(test)>2):
            if(test[0]==param and test[1]=='='):
                return float(test[2])
    return 0

def Parse(infile,param):
    inny=open(infile,'r')
    for line in inny.readlines():
        test=line.split()
        if(len(test)>2):
            if(test[0]==param and test[1]=='='):
                return (test[2])
    return 0

def findmax(array,col):
    test=float(array[0][col])
    imax=0
    for i in range(len(array)):
        if(float(array[i][col])>test):
            test=float(array[i][col])
            imax=i
    return imax

def readpeaklist(infile,met_ppm,ile_ppm):
    peak=[]
    peakfile=open(infile,'r')
    cnt=0
    for line in peakfile.readlines():
        cnt=cnt+1
        if(cnt>2):
            linetosave=string.split(line)
            peak.append(linetosave)
    peakfile.close()

    #first, work out which column is carbon and which is proton
    max1=float(peak[findmax(peak,1)][1])
    if(max1<12):
        print('column 1 is proton, column 2 is carbon')
    else:
        print('column 1 is carbon, column 2 is proton')
        peak_tmp=[]
        for i in range(len(peak)):
            peak_tmp.append((peak[i][0],peak[i][2],peak[i][1]))
        peak=peak_tmp
        
    
    
    peak_new=peak
    peak_new_full=[]
    cntILE=0
    cntMET=0
    cntLV=0
    for i in range(len(peak_new)):
        if(float(peak_new[i][1])>met_ppm):
            cntMET+=1
            line=peak_new[i][0],peak_new[i][1],peak_new[i][2],'M'
        elif(float(peak_new[i][2])<ile_ppm):
            cntILE+=1
            line=peak_new[i][0],peak_new[i][1],peak_new[i][2],'I'
        else:
            cntLV+=1
            line=peak_new[i][0],peak_new[i][1],peak_new[i][2],'LV'
        peak_new_full.append(line)

    print()
    print('Isoleucine count from chemical shifts:     ',cntILE)
    print('Leucine valine count from chemical shifts: ',cntLV)
    print('Methionine count from chemical shifts: ',cntMET)
    print()

    return peak_new_full


def peak_fold2(peak,d_min,d_max,dd):
    peak_new=[]
    for i in range(len(peak)):
        if(numpy.fabs(float(peak[i][2])-d_min)<1E-2 or numpy.fabs(float(peak[i][2])-d_max)<1E-2): #don't fold if on the edge of the spectrum and rounding errors
            peak_new.append((peak[i][0],peak[i][1],peak[i][2]))
        elif(float(peak[i][2])<=d_min):
            peak_new.append((peak[i][0],float(peak[i][1]),float(peak[i][2])+(d_max-d_min)+dd))
        elif(float(peak[i][2])>=d_max):
            peak_new.append((peak[i][0],float(peak[i][1]),float(peak[i][2])-(d_max-d_min)-dd))
        else:
            peak_new.append((peak[i][0],peak[i][1],peak[i][2]))
    return peak_new


class multiFrame(wx.Panel):
    """ The main frame of the application
    """
    title = '2D slices of 3D data'


    def __init__(self,parent):
    #def __init__(self,uc1min,uc1max,peak,index_data,thresh,offset,conn_data,spectrumfile):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        self.READ=0  #zero if data not read in, 1 if not

        self.parent=parent

        self.lc=wx.ListCtrl(self,-1,style=wx.LC_REPORT|wx.SUNKEN_BORDER)
        self.lc.InsertColumn(0, 'path')
        self.lc.InsertColumn(1, 'shize')
        self.lc.SetColumnWidth(0, 140)

        items='hnco','hnca','hncoca','hnco'
        cnt=0
        for no in items:
            self.lc.InsertStringItem(cnt,str(no))
            self.lc.SetStringItem(cnt,0,str(no))

        self.buttonRead =  wx.Button(self, -1, 'Read')

        self.vbox = wx.BoxSizer(wx.VERTICAL)     
        self.vbox.Add(self.lc, 1, wx.EXPAND)
        self.vbox.Add(self.buttonRead,1)

        #vbox=wx.BoxSizer(wx.VERTICAL)
        #vbox.Add(self.Addbutton, 0, wx.ALIGN_CENTER| wx.TOP)
        #vbox.Add(self.Removebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        #vbox.Add(self.Clearbutton, 0, wx.ALIGN_CENTER| wx.TOP)
        #vbox.Add(self.Closebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        #vbox.Add(self.Savebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        #vbox.Add(self.textbox, 0, wx.ALIGN_CENTER| wx.TOP)
        #hbox.Add(vbox)
        self.SetSizer(self.vbox)
        #vbox.Fit(self)
        #self.Centre()
        #self.Show(True)

        self.buttonRead.Bind(wx.EVT_BUTTON, self.OnButtonRead)

        sig2=ParseFlt('deconPar','sig2')
        sig3=ParseFlt('deconPar','sig3')
        indir=Parse('deconPar','indir')
        thresh=ParseFlt('deconPar','thresh')
        fac=ParseFlt('deconPar','fac')
        squash=ParseFlt('deconPar','squash')
        infile=Parse('deconPar','infile')
        peakfile=Parse('deconPar','peakfile')
        self.conn_data=[]        

        #print 'Initial peak width dimension 1: ',sig1
        print('Initial peak width dimension 2: ',sig2)
        print('Initial peak width dimension 3: ',sig3)
        print('SignalToNoise threshold:        ',thresh)


    def GetRange(self,infile):
        dic,data=ng.pipe.read(infile)
        uc0 = ng.pipe.make_uc(dic,data,dim=0)
        uc1 = ng.pipe.make_uc(dic,data,dim=1)
        uc2 = ng.pipe.make_uc(dic,data,dim=2)

        Size=data.shape
        uc0max=uc0.ppm(0)
        uc0min=uc0.ppm(Size[0]-1)
        uc1max=uc1.ppm(0)
        uc1min=uc1.ppm(Size[1]-1)
        uc2max=uc2.ppm(0)
        uc2min=uc2.ppm(Size[2]-1)
        
        print("Spectrum dimensions (pts): ",Size)   #print the spectral dimensions
        #print "Labels: ",self.labb
        print("dimension 0 limits (ppm): ", uc0min, uc0max)  #carbon 
        print("dimension 1 limits (ppm): ", uc1min, uc1max)  #direct 
        print("dimension 2 limits (ppm): ", uc2min, uc2max)  #direct 
        return uc0min,uc0max,uc1min,uc1max,numpy.fabs(uc1.ppm(0)-uc1.ppm(1))


    def OnButtonRead(self,event):
        #indir=self.dirBox.GetValue()
        #infile=self.infileBox.GetValue()
        #peakfile=self.peakBox.GetValue()

        self.a0,self.a1,self.a2,self.a3,self.aa=self.GetRange('hncoca/testJig.ft3')

        self.peak=readpeaklist('hnco'+'/'+'test.list',1.7,17.0)
        self.peak=peak_fold2(self.peak,self.a2,self.a3,self.aa)

        #uc1min,uc1max,peak,noiseVal=self.makeinp(indir,infile,peakfile)
        #self.peak=peak
        #self.uc1min=uc1min
        #self.uc1max=uc1max
        #self.noiseVal=noiseVal
        #self.spectrumfile=indir+'/'+infile
        
        if(self.READ==0):
            #print self.cb_grid.IsChecked()
            #self.parent.AddTabTwo(True,self)
            self.parent.AddTabThree(True,self)
            self.parent.AddTabFour(True,self)
            self.READ=1
