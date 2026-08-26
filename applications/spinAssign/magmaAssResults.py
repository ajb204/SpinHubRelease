#!/usr/bin/python
import numpy,sys,os,string,math,wx,glob,re
import nmrglue as ng

import wx,string,copy,math,numpy,os
import matplotlib            #import matplotlib
# matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.collections import LineCollection
import matplotlib.cm as cm
from matplotlib.figure import Figure
from wx.lib.mixins.listctrl import ColumnSorterMixin
from wx.lib.mixins.listctrl import ListCtrlAutoWidthMixin
import assign.textEdit
import matplotlib.pyplot as plt

#from magma4.decon.shiftXPostFilter import shiftXNMR
#from magma4.magma.nofit import Visualise
#import magma4.magma.analysis as analysis

from decon.shiftXPostFilter import shiftXNMR
from nofit import Visualise
import analysis as analysis


class SortedListCtrl(wx.ListCtrl, ColumnSorterMixin):
    def __init__(self, parent,dicty):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ColumnSorterMixin.__init__(self,len(dicty.keys()))
        self.itemDataMap = dicty

    def GetListCtrl(self):
        return self

class MagmaResultsMan(wx.App):
    def __init__(self,inherit):
        self.frame_ProcessFrame=MagmaResultsFrame(None,30,'MagmaResults',inherit)
        #FGA added
        self.frame_ProcessFrame.Centre(direction=wx.BOTH)
        self.frame_ProcessFrame.Show(True)
#        return Frame1(parent)

# assign ID numbers
[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
] = [wx.NewId() for _init_ctrls in range(4)]


#FGA added
class AutoWidthListCtrl(wx.ListCtrl, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        # wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT,size=(650,-1))
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ListCtrlAutoWidthMixin.__init__(self)

class MagmaResultsFrame(wx.Frame):

    def __init__(self,parent,id,title,inherit):
        #wx.Panel.__init__(self, parent=parent)
        self.parent=inherit
        self.WXV=int(wx.__version__.split('.')[0])

        #self.StyleSetForeground(wx.stc.STC_STYLE_DEFAULT,wx.Colour(230, 230, 250))
        #FGA changed
        #wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
        #      pos=wx.Point(258, 184), size=wx.Size(800, 480),
        #      style=wx.DEFAULT_FRAME_STYLE, title=u'MAGMA results ...')
        #self.SetClientSize(wx.Size(900, 280))
        monitorWidth, monitorHeight = wx.GetDisplaySize()
        #wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
        #      pos=wx.DefaultPosition, size=(monitorWidth, monitorHeight),
        #      style=wx.DEFAULT_FRAME_STYLE, title=u'MAGMA results')
        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=parent,
                          pos=wx.DefaultPosition, size=(monitorWidth, monitorHeight),
                          title=u'MAGMA results')

        #self.SetBackgroundColour('WHITE')
        #self.panel=wx.Panel(self,-1)

        ########

        self.modes=wx.ComboBox(self, -1, size=(150, -1), style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.OnRefresh, self.modes)
        self.GetModes(curr='y') #set value based on current progress #check latest

        self.subgraph=wx.ComboBox(self, -1, size=(80, -1), style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.OnTickPlotFilt, self.subgraph)
        self.SetSubgraph() #set combo box based on current selection

        ########


        #self.source = wx.ListCtrl(self, -1, style = wx.LC_REPORT,size=(650,300))
        #self.source = wx.SortedListCtrl(self, -1, style = wx.LC_REPORT,size=(650,300))
        #self.lc=SortedListCtrl(panel,self.corrDict)

        self.source = AutoWidthListCtrl(self)
        self.source.SetMinSize((650,300))

        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnDoubleClick, self.source)

        #FGA commented- don't think these do anything
        self.statusLbl1  = wx.StaticText(self, label="")
        self.statusLbl2  = wx.StaticText(self, label="")


        self.confLbl  = wx.StaticText(self, label="")
        self.solnLbl  = wx.StaticText(self, label="")
        self.nmreLbl  = wx.StaticText(self, label="")
        self.pdbeLbl  = wx.StaticText(self, label="")

        lblList = ['Progress','Optimise','Distance', 'Intensity','Talos']
        self.rbox = wx.RadioBox(self,label = '',choices = lblList ,majorDimension = 0, style = wx.RA_SPECIFY_ROWS)
        self.Bind(wx.EVT_RADIOBOX,self.OnTickPlot,self.rbox)

        self.buttonClose = wx.Button(self, label="Close")
        self.Bind(wx.EVT_BUTTON, self.OnButtonClose, self.buttonClose)

        self.buttonLog = wx.Button(self, label="Log")
        self.Bind(wx.EVT_BUTTON, self.OnButtonLog, self.buttonLog)

        self.buttonSave = wx.Button(self, label="Save")
        self.Bind(wx.EVT_BUTTON, self.OnButtonSave, self.buttonSave)

        self.buttonShiftX = wx.Button(self, label="ShiftX")
        self.Bind(wx.EVT_BUTTON, self.onShiftXBtn, self.buttonShiftX)

        self.buttonRefresh = wx.Button(self, label="Refresh")
        self.Bind(wx.EVT_BUTTON, self.OnRefresh, self.buttonRefresh)

        sheetList = ['Assignments','Filters','NOEs',]
        self.sheetbox = wx.RadioBox(self,label = '',choices = sheetList ,majorDimension = 0, style = wx.RA_SPECIFY_ROWS)
        self.Bind(wx.EVT_RADIOBOX,self.OnTickFilt,self.sheetbox)


        self.buttonShow = wx.Button(self, label="ShowData")
        self.Bind(wx.EVT_BUTTON, self.OnButtonShow, self.buttonShow)

        #self.PopulateList()


        self.sizerG= wx.BoxSizer(wx.VERTICAL)
        self.sizerG.Add(self.modes)
        self.sizerG.Add(self.confLbl)
        self.sizerG.Add(self.solnLbl)
        self.sizerG.Add(self.nmreLbl)
        self.sizerG.Add(self.pdbeLbl)


        self.sizerP= wx.BoxSizer(wx.HORIZONTAL)
        self.sizerPa= wx.BoxSizer(wx.VERTICAL)
        self.sizerPa.Add(self.rbox)
        self.sizerPa.Add(self.buttonClose)
        self.sizerPa.Add(self.buttonLog)
        self.sizerPa.Add(self.buttonSave)
        self.sizerPb= wx.BoxSizer(wx.VERTICAL)
        self.sizerPb.Add(self.sheetbox)
        self.sizerPb.Add(self.buttonShiftX)
        self.sizerPb.Add(self.buttonRefresh)
        self.sizerPb.Add(self.buttonShow)

        self.sizerPb.Add(self.subgraph)


        self.sizerP.Add(self.sizerPa)
        self.sizerP.Add(self.sizerPb)



        self.statusList=[]
        self.statusBox=wx.BoxSizer(wx.VERTICAL)
        #FGA changed to add r6 mode
        for i in range(10):
            self.statusList.append(wx.StaticText(self, label=""))
            self.statusBox.Add(self.statusList[i])






        self.sumLbl = wx.StaticBox(self, -1, 'Summary:', size=(240, 140))
        self.sumSizer = wx.StaticBoxSizer(self.sumLbl, wx.VERTICAL)
        self.borderS = wx.BoxSizer()
        self.borderS.Add(self.sizerG, 1, wx.ALL | wx.EXPAND, 7)
        self.sumSizer.Add(self.borderS,wx.EXPAND)


        self.pltLbl = wx.StaticBox(self, -1, 'Options:', size=(240, 140))
        self.pltSizer = wx.StaticBoxSizer(self.pltLbl, wx.VERTICAL)
        self.borderP = wx.BoxSizer()
        self.borderP.Add(self.sizerP, 1, wx.ALL | wx.EXPAND, 7)
        self.pltSizer.Add(self.borderP,wx.EXPAND)


        self.statusLbl = wx.StaticBox(self, -1, 'Status:', size=(240, 140))
        self.statusSizer = wx.StaticBoxSizer(self.statusLbl, wx.VERTICAL)
        self.borderT = wx.BoxSizer()
        self.borderT.Add(self.statusBox, 1, wx.ALL | wx.EXPAND, 7)
        self.statusSizer.Add(self.borderT,wx.ALL|wx.EXPAND)


        self.resLbl = wx.StaticBox(self, -1, 'Results:', size=(240, 140))
        self.resSizer = wx.StaticBoxSizer(self.resLbl, wx.VERTICAL)
        self.borderR = wx.BoxSizer()
        self.borderR.Add(self.source, 1, wx.ALL | wx.EXPAND, 7)
        self.resSizer.Add(self.borderR,wx.EXPAND)


        #self.assLbl = wx.StaticBox(self, -1, 'Assignments:', size=(240, 140))
        #self.assSizer = wx.StaticBoxSizer(self.assLbl, wx.VERTICAL)
        #self.borderA = wx.BoxSizer()
        #self.borderA.Add(self.sizerA, 1, wx.ALL | wx.EXPAND, 7)
        #self.assSizer.Add(self.borderA)

        self.splitSizer = wx.GridBagSizer(1, 3)

        self.splitSizer.Add(self.sumSizer,(0,0))
        self.splitSizer.Add(self.pltSizer,(0,1))
        self.splitSizer.Add(self.statusSizer,(0,2))


        self.longSizer=wx.GridBagSizer(2, 1)

        self.longSizer.Add(self.splitSizer,(0,0))
        #FGA changed
        #self.longSizer.Add(self.resSizer,(1,0))
        self.longSizer.Add(self.resSizer,(1,0), flag=wx.EXPAND|wx.ALL)

        self.fullSizer=wx.BoxSizer(wx.HORIZONTAL)
        self.fullSizer.Add(self.longSizer)


        self.create_main_panel()
        self.draw_figure()
        self.canvas.draw()

        self.SetSizerAndFit(self.fullSizer)

        #FGA added
        #if self.ref=='subgraphMode':
        #    self.rbox.Disable()

        self.Show(True)
        self.Fit()
        #self.Freeze()
        self.Layout()
        #self.Thaw()

        #topPanel.SetSizer(sizer)
        #for child in self.panel2.GetChildren():
        #    child.Destroy()

        #print(if calc is complete or still running.)


    #if we have intensity plot and NOE tab selected,
    #then a click on the intensity graph selects an NOE.
    def on_pick(self, event):
        seleRbox=self.rbox.GetSelection()
        if(seleRbox==3):
            print('Intensity')
            print(event.xdata,event.ydata)
            x_min,x_max=self.ax.get_xlim()
            y_min,y_max=self.ax.get_ylim()
            xdist=x_max-x_min
            ydist=y_max-y_min
            rad2=((self.inty[:,0]-event.xdata)/xdist)**2.+((self.inty[:,2]-event.ydata)/ydist)**2.

            argmin=numpy.argmin(rad2)
            print('argmin:',argmin)
            print(self.sheetbox.GetStringSelection())
            if(self.sheetbox.GetStringSelection()=='NOEs'):

                sele=self.source.GetFirstSelected()
                self.source.Select(sele,on=0)
                self.source.Select(argmin,on=1)
                self.source.EnsureVisible(argmin)
                self.OnTickPlot(True)
                #self.subgraph.SetSelection(argmin)

            #self.xv = [self.source.GetItem(row, 7).GetText() for row in xrange(count)][sele]
            #self.yv = [self.source.GetItem(row, 3).GetText() for row in xrange(count)][sele]
                    #print(xv,yv)
            #self.ax.scatter(float(self.xv),float(self.yv),color='r',s=200)


    def onShiftXBtn(self,event):
        if 'shiftXFile' not in self.parent.pars.keys():
        # if not self.parent.pars['shiftXFile']:
            print('Could not compare with calculated NMR shifts as the file is not there.')
            self.sheetbox.SetStringSelection('Assignments')
            self.OnTickFilt(True)
            return
        if 'ileShiftFile' not in self.parent.pars.keys():
        # if not self.parent.pars['ileShiftFile']:
            print('Could not compare with calculated NMR shifts as the experimental NMR shift file is not there.')
            self.sheetbox.SetStringSelection('Assignments')
            self.OnTickFilt(True)
            return

        self.doShiftX()
        self.PopulateList(shift='y')

        """
        #populate list
        resultConf={}
        conf=0
        for key,vals in resDict.items():
            if(len(vals)==1):
                conf+=1
                resultConf[key]=vals
                del(resDict[key])
        cnt=0
        for i in range(self.colVal+1):
            if(os.path.exists(self.progress[i][1]+'/confident.res')): #is this mode done?
                new={}
                inny=open(self.progress[i][1]+'/confident.res')
                for line in inny.readlines():
                    test=line.split(':')
                    key=test[0].split()[0]
                    ass=test[1].split()[0]
                    if(key in resultConf.keys()):
                        if(ass==resultConf[key][0]):
                            new[key]=ass
                            del(resultConf[key])
                inny.close()

                for key,vals in new.items():
                    cnt+=1

                    num_items = self.source.GetItemCount()
                    self.source.InsertStringItem(num_items,str(cnt))
                    self.source.SetStringItem(num_items,0,str(cnt))
                    self.source.SetStringItem(num_items,1,key)
                    self.source.SetStringItem(num_items,2,self.progress[i][0])

                    #num_items = self.source.GetItemCount()
                    ##FGA changed- depreciated functions
                    #self.source.InsertItem(str(cnt))
                    #self.source.SetItem(0,str(cnt))
                    #self.source.SetItem(1,key)
                    #self.source.SetItem(2,self.progress[i][0])


                    stry=vals
                    self.source.SetStringItem(num_items,3,stry)

                    if(self.progress[i][0]=='subgraphMode' or self.progress[i][0]=='polishMode' or self.progress[i][0]=='nudgeMode' or self.progress[i][0]=='finalMode'):

                        color = (0,int(255), 0)
                        self.source.SetItemBackgroundColour(num_items,color)

        for key in sorted(resultConf,key=lambda k:len(resultConf[k])):
            val = resultConf[key][0]
            cnt+=1
            num_items = self.source.GetItemCount()
            self.source.InsertItem(num_items,str(cnt))
            self.source.SetItem(num_items,0,str(cnt))
            self.source.SetItem(num_items,1,key)
            self.source.SetItem(num_items,2,'compareShiftsFilter')
            self.source.SetItem(num_items,3,val)

        for key in sorted(resDict,key=lambda k:len(resDict[k])):
            vals=resDict[key]
            cnt+=1
            num_items = self.source.GetItemCount()
            self.source.InsertItem(num_items,str(cnt))
            self.source.SetItem(num_items,0,str(cnt))
            self.source.SetItem(num_items,1,key)
            stry=''
            for val in vals:
                stry+=val+' '
            self.source.SetItem(num_items,3,stry)

        print('Comparison with calculated shifts done.')
        """

    def doShiftX(self):
        print('shiftXPostMode: Ruling out options from combined results if we can.')

        inst = self.parent.inst
        #Magma('input.magma',run='n')
        prior = str(self.colVal+1)
        if os.path.exists(inst.P.outdir+'/'+prior):
            shiftX = shiftXNMR(self.parent.pars['shiftXFile'].strip(), self.parent.pars['ileShiftFile'].strip(),
                                                self.parent.pars['chains'], self.parent.pars['residues'],
                                                inst.P.outdir+'/'+prior+'/combinedResults.res')
        else:
            shiftX = shiftXNMR(self.parent.pars['shiftXFile'].strip(), self.parent.pars['ileShiftFile'].strip(),
                                                self.parent.pars['chains'], self.parent.pars['residues'],
                                                inst.P.outdir+'/combinedResults.res')
        self.shiftDict=shiftX.Parse() #return assignment dictionary







    def SetSubgraph(self):
        #self.subgraph=wx.ComboBox(self, -1, size=(80, -1), choices=listy, style=wx.CB_READONLY)
        #self.Bind(wx.EVT_COMBOBOX, self.OnTickPlot, self.subgraph)
        listy=[]
        listy.append('all')
        if(self.ref=='subgraphMode' or self.ref=='polishMode'):
            for key,vals in self.parent.inst.subgraphRef.items():
                listy.append(str(key+1))
        if(self.ref=='nudgeMode'):
            self.grps=self.DoNudge()
            for i in range(len(self.grps)):
                listy.append(str(i))

        self.subgraph.SetItems(listy)

        #FGA changed
        #if self.ref=='nudgeMode':
        #    self.subgraph.SetSelection(1)
        #else:
        #    self.subgraph.SetSelection(0)
        self.subgraph.SetSelection(0)


    def ParseMagma(self):
        inny=open(self.parent.inst.P.outdir+'/input.magma')
        check=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                if(test[0]=='subgraphMode' and test[1]=="on"):
                    check.append(0)
                if(test[0]=='polishMode' and test[1]=="on"):
                    check.append(1)
                if(test[0]=='nudgeMode' and test[1]=="on"):
                    check.append(2)
                if(test[0]=='finalMode' and test[1]=="on"):
                    check.append(3)
                if(test[0]=='longMode' and test[1]=="on"):
                    check.append(4)
        if(self.WXV==4):
            self.parent.protoBox.SetCheckedItems(check)
        else:
            self.parent.protoBox.SetChecked(check)


    def GetModes(self,curr='n'):
        if(curr=='n'):
            self.colVal=self.modes.GetCurrentSelection()

        self.ParseMagma()
        self.progress=[]

        cnt=0
        if(self.parent.protoBox.IsChecked(0)):
            cnt+=1
            self.progress.append(('subgraphMode',self.parent.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.parent.protoBox.IsChecked(1)):
            cnt+=1
            self.progress.append(('polishMode',self.parent.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.parent.protoBox.IsChecked(2)):
            cnt+=1
            self.progress.append(('nudgeMode',self.parent.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.parent.protoBox.IsChecked(3) or self.parent.protoBox.IsChecked(0)==False):
            cnt+=1
            self.progress.append(('finalMode',self.parent.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))
        if(self.parent.protoBox.IsChecked(4)):
            cnt+=1
            self.progress.append(('longMode',self.parent.inst.P.outdir+'/'+str(cnt),cnt,'incomplete'))


        self.progress=numpy.array(self.progress)
        self.progress[-1][1]=self.parent.inst.P.outdir #adjust last folder

        cnt=0
        for i in range(len(self.progress)-1):
            if(os.path.exists(self.progress[i][1])):
                self.progress[i,3]='complete'
                cnt+=1
        if(cnt==len(self.progress)-1):
            if(os.path.exists(self.progress[-1,1]+'/combinedResults.res')):
                self.progress[-1,3]='complete'

        #print(self.progress)
        #print(len(self.progress))
        tmp=[]
        tag=0
        for i in range(len(self.progress)):
            if(self.progress[i,3]=='incomplete'):
                tmp.append(self.progress[i,:])
                break
            else:
                tmp.append(self.progress[i,:])
        self.progress=numpy.array(tmp)
        print('after',len(self.progress))
        print(self.progress)
        self.modes.SetItems(self.progress[:,0])

        if(curr=='y'):
            self.colVal=len(self.progress)-1
            while(1==1):
                if(self.colVal>0):
                    if(self.progress[self.colVal-1,3]=='incomplete'):
                        self.colVal-=1
                    else:
                        break
                else:
                    break

        self.modes.SetSelection(self.colVal)
        self.ref = self.progress[self.colVal][0]
        self.testdir=self.progress[self.colVal][1]
        if(os.path.exists(self.testdir)==0):
            self.testdir=self.parent.inst.P.outdir
            print('Calculation not finished.')
            print('Current mode is not saved')


    def ParseLogAll(self,var,full='n',split='y'):
        inny=open(self.testdir+'/log')
        listy=[]
        tag=0
        for line in inny.readlines():
            test=line.split(var)
            if(len(test)>1):
                tag=1
                if(full=='n'):
                    if(split=='y'):
                        return line.split(var)[1].split()
                    else:
                        return line.split(var)[1]
                else:
                    if(split=='y'):
                        listy.append(line.split(var)[1].split())
                    else:
                        listy.append(line.split(var)[1])
        if(tag==0):
            print('Could not find',var,'in')
            return 0
        else:
            return listy

    def DoNudge(self):
        inny=open(self.testdir+'/log')
        iny=0
        grp=0
        tbl=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)!=0):
                if(test[0]=='INTERSECT:'):
                    iny=1
                if(test[0]=='Groupings:'):
                    grp=1
                if(iny==1 and grp==0 and test[0]!='INTERSECT:'):
                    tbl.append(test)
        inny.close()
        #st=''
        #for ii in range(len(self.inst.subgraphRef.keys())):
        #    st+='c'
        #outlat.write('\\begin{tabular}{c|%s}\n' % (st))
        #for ii in range(len(self.inst.subgraphRef.keys())):
        #    outlat.write('& %i' % (ii+1))
        #outlat.write('\\\\ \n')
        #outlat.write('\\hline\n')
        #for ii in range(len(self.inst.subgraphRef.keys())):
        #    outlat.write(' %i' % (ii+1))
        #    for jj in range(len(self.inst.subgraphRef.keys())):
        #        outlat.write('& %s' % (tbl[ii][jj+1].split(':')[1]))
        #    outlat.write('\\\\ \n')
        #outlat.write('\\end{tabular}\n')
        #outlat.write('\n\n')
        #outlat.write('Grouping subgraphs with overlapped assignments:\n\n ')

        groups=self.ParseLogAll("Groupings:",split='n')
        test= re.findall(r'\d+',groups)

        grps=[]
        cnt=0
        proggy=''
        for ii,t in enumerate(test):
            proggy+=str(int(t)+1)+' '
            if(ii%self.parent.inst.P.nudgy==self.parent.inst.P.nudgy-1):
                grps.append((cnt,proggy))
                proggy=''
                cnt+=1
        return grps


    def WhereAmI(self): #figure out where in the calc we are.
        #self.testdir=self.progress[self.colVal][1]
        print('working out calculation stage...')
        self.status=0

        #if save directory does not exist, calc is incomplete
        #if combinedResults.res does not exist, calc is incomplete
        for i in range(5):
            self.statusList[i].SetLabel("")

        self.statusList[0].SetLabel('mode '+self.progress[self.colVal,3])

        stt=2
        #check for existance of progress files
        if(self.ref=='subgraphMode' or self.ref=='polishMode'):
            wait=[]
            win=[]
            for ref in self.parent.inst.subgraphRef.keys():
                vf2file=self.testdir+'/vf2_%i.txt' % (int(ref)+1)
                mcsfile=self.testdir+'/mces_%i.txt' % (int(ref)+1)
                tag=0
                if(os.path.exists(vf2file)):
                    win.append(str(int(ref)+1)+'vf2' )
                    tag=1
                if(os.path.exists(mcsfile)):
                    win.append(str(int(ref)+1)+'mcs' )
                    tag=1
                if(tag==0):
                    wait.append(str(int(ref)+1))
                #else:

            strg=''
            strg='Finished: '
            for i in range(len(win)):
                strg+=str(win[i])+' '
                if(len(strg)>25):
                    self.statusList[stt].SetLabel(strg)
                    strg=''
                    stt+=1
            self.statusList[stt].SetLabel(strg);stt+=1

            if(len(wait)>0):
                strg=''
                strg='Remain: '
                for i in range(len(wait)):
                    strg+=str(wait[i])+' '
                self.statusList[stt].SetLabel(strg)

            if(len(wait)>0):
                print('currently running ',wait[0])
                curr=' on '+str(wait[0])
                if(len(wait)>1):

                    self.status=1
                    print('to go...')
                    for i,wa in enumerate(wait):
                        if(i!=0):
                            print(wa)
            else:
                curr=''
        elif(self.ref=='nudgeMode'):
            self.grps=self.DoNudge()
            stt=2
            strg='groups: '
            for i in range(len(self.grps)):
                strg+='('+self.grps[i][1]+') '
                if(len(strg)>25):
                    #print(stt)
                    self.statusList[stt].SetLabel(strg)
                    strg=''
                    stt+=1
            self.statusList[stt].SetLabel(strg)
        else:
            vf2file=self.testdir+'/vf2.txt'
            mcsfile=self.testdir+'/mces.txt'
            tag=0
            if(os.path.exists(vf2file)):
                tag=1
            if(os.path.exists(mcsfile)):
                tag=1
            #if(tag==0):
            curr=''

        self.GetLastFile()
        #print('fil:',self.fil)
        if(self.fil==''):
            return

        if(len(self.fil.split('community'))>1):
            self.status=2
            self.statusList[1].SetLabel('running communityMode %s' % curr)
        if(len(self.fil.split('_new'))>1):
            self.status=3
            self.statusList[1].SetLabel('running splitMode %s' % curr)
            return
        if(len(self.fil.split('vf2'))>1):
            self.status=4
            self.statusList[1].SetLabel('running vf2')
            return
        if(len(self.fil.split('.G'))>1):
            self.status=5
            self.statusList[1].SetLabel('running mces %s' % curr)
            return
        if(len(self.fil.split('conv'))>1):
            self.status=6
            self.statusList[1].SetLabel('optimising mces %s' % curr)
            return

    def GetLastFile(self):
        list_of_files = glob.glob(self.testdir+'/*') # * means all if need specific format then *.csv
        adjusted= sorted(list_of_files,key=os.path.getctime)
        #print(adjusted)
        #latest_file = max(list_of_files, key=os.path.getctime)
        #print(latest_file)
        cnt=0
        while(1==1):
            ii=len(adjusted)-1-cnt
            test=adjusted[ii].split('/')[-1]
            #print(test)
            if(test!='mces_traj.png' and test!='mcesplt.gp' and test!='log'):
                self.fil=test
                return
            cnt+=1
            if(cnt==len(adjusted)):
                self.fil=''
                return

    def AtoI(self,val):
        for i in range(len(self.parent.parent.tabOne.peak)):
            if(val==self.parent.parent.tabOne.peak[i].name):
                return i

    def OnDoubleClick(self,event):
        self.OnTickPlot(True)
        print(self.sheetbox.GetStringSelection())
        if(self.sheetbox.GetStringSelection()=='NOEs'):
                sele=self.source.GetFirstSelected()
                print(sele)
                count = self.source.GetItemCount()
                col1 = [self.source.GetItem(row, 1).GetText() for row in xrange(count)][sele]
                col2 = [self.source.GetItem(row, 2).GetText() for row in xrange(count)][sele]
                print(col1,col2)

                self.parent

                self.parent.parent.tabTwo.ComboBox1.SetSelection(self.parent.parent.tabTwo.AtoI(col2))
                self.parent.parent.tabTwo.ComboBox2.SetSelection(self.parent.parent.tabTwo.AtoI(col1))

                self.parent.parent.tabTwo.draw_figure()

                #self.parent.parent.tabFour.ComboBox1.SetSelection(self.AtoI(p1))
                #self.parent.parent.tabFour.ComboBox2.SetSelection(self.AtoI(p2))
                #self.parent.parent.tabFour.draw_figure()



    def OnButtonClose(self,event):
        self.Close()

    def OnButtonLog(self,event):
        textEdit.MyFrame(self.testdir+'/log')

    def OnButtonSave(self,event):
        file_choices='*'
        dlg = wx.FileDialog(
            self,
            message="Save session...",
            defaultDir=os.getcwd(),
            defaultFile='results.res',
            wildcard=file_choices,
            style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()

            print('Writing results to file: ',path)
            outfileFile=path

            #if(os.path.exists(outfileFile)==0):
            #    outy=open(outfileFile,'w');outy.close()
            os.system('cp '+self.parent.inst.P.outdir+'/combinedResults.res '+outfileFile)

            self.UpdateScores()
            self.parent.parent.tabOne.molecule.assemble(outfileFile+'.peaks',self.result_dict)



            #inst = self.parent.inst
            #Magma('input.magma',run='n')
            #prior = str(self.colVal+1)
            #if os.path.exists(inst.P.outdir+'/'+prior):

            #self.canvas.print_figure(path, dpi=self.dpi)
            #self.parent.parent.flash_status_message("Saved %s" % path)


    def OnButtonShow(self,event):
        test=Visualise(self.parent.inst)

        tig=0
        if(self.ref=='subgraphMode' or self.ref=='polishMode'):
            sele=self.subgraph.GetStringSelection()
            if(sele=='all'):
                tig=1
        else:
            tig=1

        if(tig==1):
            test.monte_min()
            test.pymol_gen_noe()
        else:
            test.monte_min(sele=int(sele))
            test.pymol_gen_noe()

        #from os import popen
        self.parent.PymolExec('pyscript_noe.py')
        #self.calcy=popen('pymol report/pyscript_noe.py')


    def OnTickFilt(self,event):
        self.WhereAmI() #figure out where in the calc we are.
        #selebox=self.sheetbox.GetSelection()
        boxVal=self.sheetbox.GetStringSelection()
        #print(selebox,self.sheetbox.GetStringSelection())
        self.source.ClearAll()
        if(boxVal=='Assignments'):
            self.source.InsertColumn(0, 'id', width = 50,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(1, 'NMR', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(2, 'mode', width = 100,format=wx.LIST_FORMAT_CENTRE)

            self.source.InsertColumn(3, 'PDB',width=1000)
            self.PopulateList()
            return
        elif(boxVal=='Filters'):
            self.source.InsertColumn(0, 'id', width = 50,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(1, 'NMR', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(2, 'mode', width = 100,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(3, 'PDB', width = 10000)
            self.PopulateFilt()
            return
        elif(boxVal=='NOEs'):
            self.source.InsertColumn(0, 'id', width = 50,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(1, 'A', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(2, 'B', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(3, 'intensity', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(4, '#', width = 50,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(5, 'tp', width = 50,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(6, 'hits', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(7, 'dist', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(8, 'stdev', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.PopulateNOE()
            return
        elif boxVal == 'CompareShifts':
            self.source.InsertColumn(0, 'id', width = 50,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(1, 'NMR', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(2, 'mode', width = 100,format=wx.LIST_FORMAT_CENTRE)
            #FGA changed
            # self.source.InsertColumn(3, 'PDB', width = 10)
            self.source.InsertColumn(3, 'PDB')
            self.onShiftXBtn()


    def OnRefresh(self,event):
        print('getting modes')
        self.GetModes()
        print('setting subgraph')
        self.SetSubgraph() #reset subgraph options
        self.parent.UpdateReport()
        self.sheetbox.SetSelection(0)
        self.OnTickPlotFilt(True)
        #FGA added
        #self.Layout()

    def OnTickPlotFilt(self,event):
        self.UpdateScores()
        self.WhereAmI()
        self.OnTickPlot(True)
        self.OnTickFilt(True)

    def OnTickPlot(self,event):
        self.draw_figure()
        #FGA added
        #if self.ref=='subgraphMode':
        #    self.rbox.Disable()
        #else:
        #    self.rbox.Enable()


    def create_main_panel(self):

        self.fig=Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.mpl_connect('button_press_event', self.on_pick)

        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)


        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        #self.vbox.AddSpacer(10)
        #self.SetSizer(self.vbox)
        self.vbox.Fit(self)
        self.fullSizer.Add(self.vbox)

    #read in G2hist file for distance and intensity plots
    #combine results for subgraph/polishmodes
    def GetDistHist(self):
        if(self.ref=='subgraphMode' or self.ref=='polishMode' or self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()
            if(sele=='all'):
                #print('Reading all subgraphs')
                for ii in range(len(self.parent.inst.subgraphRef.keys())):
                    testfile=self.testdir+'/G2hist_%i.hist' % (ii+1)
                    if(os.path.exists(testfile)):
                        histA,badNOEA,intyA=self.GetDistHistCore(testfile)
                        try:
                            hist[:,2]+=histA[:,2]

                            for n in intyA:
                                inty.append(n)
                            for n in badNOEA:
                                    badNOE.append(n)
                        except:
                            hist=copy.deepcopy(histA)
                            badNOE=copy.deepcopy(badNOEA)
                            inty=copy.deepcopy(intyA)
                        #print(testfile,len(inty),len(intyA))
            else:
                #for ii in range(len(self.parent.inst.subgraphRef.keys())):
                    print('Reading subgraph',sele)
                    testfile=self.testdir+'/G2hist_%s.hist' % (sele)
                    if(os.path.exists(testfile)):
                        histA,badNOEA,intyA=self.GetDistHistCore(testfile)
                        try:
                            hist[:,2]+=histA[:,2]

                            for n in intyA:
                                inty.append(n)
                            for n in badNOEA:
                                    badNOE.append(n)
                        except:
                            hist=copy.deepcopy(histA)
                            badNOE=copy.deepcopy(badNOEA)
                            inty=copy.deepcopy(intyA)
                        #print(testfile,len(inty),len(intyA))
        else:
            if(os.path.exists(self.testdir+'/G2hist.hist')):
                hist,badNOE,inty=self.GetDistHistCore(self.testdir+'/G2hist.hist')
            else:
                print('Cannot find',self.testdir+'/G2hist.hist')
        try:
            for i in range(len(hist)):
                if(hist[i,1]!=0):
                    hist[i,3]=hist[i,2]/hist[i,1]
                else:
                    hist[i,3]=0.
            return hist,badNOE,numpy.array(inty)
        except:
            return 0,0,0


    def GetDistHistCore(self,infile):
        tmp=[]
        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            tmp.append(test)
        cnt=0
        badNOE=[]
        hist=[]
        inty=[]
        for i in range(len(tmp)):
            try: #check if this line, and the next one are empty
                if(len(tmp[i])==0  and len(tmp[i+1])):
                    cnt+=1
            except:
                pass

            if(cnt==1):
                try:
                    if(int(tmp[i][5])==0):
                        badNOE.append((tmp[i][1],tmp[i][2]))
                except:
                    pass

            if(cnt==2 and len(tmp[i])!=0):
                dist=float(tmp[i][0])
                pdb=float(tmp[i][1])
                soln=float(tmp[i][2])
                ratio=float(tmp[i][3])
                hist.append((dist,pdb,soln,ratio))

            if(cnt==1 and len(tmp[i])!=0):
                dist=float(tmp[i][9])
                disterr=float(tmp[i][10])
                iV=float(tmp[i][7])
                iVerr=float(tmp[i][8])
                if(disterr=='nan'):
                    disterr=0.
                if(disterr!=disterr):
                    disterr=0.

                #print(dist,disterr,iV,iVerr)
                inty.append((dist,disterr,iV,iVerr))

        hist=numpy.array(hist)
        return hist,badNOE,inty


    def GetDistHistNOE(self,infile):
        dat=[]
        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            dat.append(test)
        inny.close()
        cnt=0
        for i,tmp in enumerate(dat):
            try: #check if this line, and the next one are empty
                if(len(tmp)==0  and len(dat[i+1])):
                    cnt+=1
            except:
                pass

            if(cnt==1 and len(tmp)>0):
                row=[]
                #row.append(int(tmp[0])) #index
                row.append(len(self.noe)) #index
                row.append(tmp[1]) #label1
                row.append(tmp[2]) #label2
                row.append(tmp[3]) #type
                row.append(int(tmp[4])) #weight
                row.append(int(tmp[5])) #number
                row.append(tmp[6]) #weighting
                row.append(tmp[7]) #intensity
                row.append(tmp[8]) #scr
                row.append('%.3f' % float(tmp[9])) #distance

                if(tmp[10]=='nan'):
                    row.append(0) #stdev
                else:
                    if(float(tmp[10])<0.1):
                        row.append(0) #stdev
                    else:
                        row.append('%.3f'% float(tmp[10])) #stdev


                self.noe.append(row)

    def GetConvFile(self,inconv1,inconv2):
        inny=open(inconv1)
        in1=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)==3):
                in1.append((float(test[0]),float(test[1]),float(test[2])))
        inny.close()

        inny=open(inconv2)
        in2=[]
        tost=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)==3):
                tost.append((float(test[0]),float(test[1]),float(test[2])))
            else:
                if(len(tost)>0):
                    in2.append(numpy.array(tost))
                tost=[]
        if(len(tost)>0):
            in2.append(numpy.array(tost))
        inny.close()
        return in1,in2


    def GetConv(self):
        analdir=self.progress[self.colVal][1]
        if(os.path.exists(analdir)==0):
            analdir=self.parent.inst.P.outdir
        #print(analdir)
        if(self.ref=='subgraphMode' or self.ref=='polishMode' or self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()
            if(sele!='all'):
                inconv1=analdir+'/conv.1_%i.out' % int(sele)
                inconv2=analdir+'/conv.2_%i.out' % int(sele)
                if(os.path.exists(inconv1)==1):
                    return self.GetConvFile(inconv1,inconv2)
        else:
                inconv1=analdir+'/conv.1.out'
                inconv2=analdir+'/conv.2.out'
                if(os.path.exists(inconv1)==1):
                    return self.GetConvFile(inconv1,inconv2)
        return 0,0

    def draw_figure(self):
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        self.ax.clear()

        seleRbox=self.rbox.GetSelection()
        #print(seleRbox,self.rbox.GetStringSelection())
        if(seleRbox==2):
            hist,badNOE,inty=self.GetDistHist()
            try:
                if(hist==0):
                    print('No G2hist file')
                    return
            except:
                pass
            #print(hist)
            #setup matplotlib


            #my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
            #plt.grid(True)
            self.ax2=self.ax.twinx()
            self.ax2.set_ylabel('Ratio',color='r')
            self.ax2.set_ylim([0.0,1.05])

            self.ax.set_xlim([3.8,20])
            self.ax2.set_xlim([3.8,20])

            l1=self.ax.plot(hist[:,0],hist[:,1],label='PDB')
            l2=self.ax.plot(hist[:,0],hist[:,2],label='Data')
            l3=self.ax2.plot(hist[:,0],hist[:,3],color='r',label='ratio = Data/PDB')

            self.ax.fill_between(hist[:,0],0,hist[:,1],facecolor='blue',alpha=0.1)
            self.ax2.fill_between(hist[:,0],0,hist[:,3],facecolor='red',alpha=0.1)
            self.ax.fill_between(hist[:,0],0,hist[:,2],facecolor='green',alpha=0.1)

            self.ax.set_xlabel('C-C Distance (A)')
            self.ax.set_ylabel('Counts')
            #ax.set_title(analdir)

            ymin,ymax=self.ax.get_ylim()
            self.ax.set_ylim([0,ymax])
            #print(ymin,ymax)

            #ax.legend((l1,l2,l3),('PDB','Data','Ratio'),'upper right')
            #ax.legend()
            #ax2.legend()

            import matplotlib.patches as mpatches
            red_patch = [mpatches.Patch(color='red', label='ratio')]
            blue_patch = [mpatches.Patch(color='blue', label='PDB')]
            green_patch = [mpatches.Patch(color='green', label='Data')]

            handles,labels=self.ax.get_legend_handles_labels()
            #handles=#handles+red_patch
            handles=blue_patch+green_patch+red_patch
            labels=labels+[u'ratio']
            self.ax.legend(handles,labels,loc=1)

            #fig.legend(handles=
            #ax.legend(('r','b','g'),("A","B","C"),loc=1)
            self.ax2.tick_params('y',colors='r')



        elif(seleRbox==3):
            hist,badNOE,inty=self.GetDistHist()
            try:
                if(hist==0):
                    print('No G2hist file')
                    return
            except:
                pass

            #print(inty)
            #print(len(inty))

            self.inty=inty
            self.ax.errorbar(inty[:,0], inty[:,2], xerr=inty[:,1], yerr=inty[:,3],fmt='o')
            yvals=numpy.ones_like(inty[:,0])*0.5


            if(self.sheetbox.GetSelection()!=2):
                self.sheetbox.SetSelection(2)
                self.OnTickFilt(True)
            if(self.sheetbox.GetStringSelection()=='NOEs'):
                #try:
                    sele=self.source.GetFirstSelected()
                    count = self.source.GetItemCount()
                    xv = [self.source.GetItem(row, 7).GetText() for row in xrange(count)][sele]
                    yv = [self.source.GetItem(row, 3).GetText() for row in xrange(count)][sele]

                    #print(xv,yv)
                    self.ax.scatter(float(xv),float(yv),color='r',s=200)

            self.ax.plot(inty[:,0],yvals)

            self.ax.set_xlabel("distance(A)")
            self.ax.set_ylabel("CrosspeakIntensity")


            xmin,xmax=self.ax.get_xlim()
            ymin,ymax=self.ax.get_ylim()

            #print(inty[:,0])
            print('xmax:',numpy.max(inty[:,0]))

            self.ax.set_xlim(0,numpy.max(inty[:,0])*1.05)
            self.ax.set_ylim(0,1.1)


        elif(seleRbox==0):
            self.OnTickFilt(True)
            try:
                prog,title=self.GetProg()
            except:
                self.fig.clear()
                print('No mces file to plot')
                return

            try:
                testmax=0
                for pro in prog:
                    test=numpy.max(pro[:,2])
                    if(test>testmax):
                        testmax=test
                testmin=testmax
                for pro in prog:
                    test=numpy.min(pro[:,2])
                    if(test<testmin):
                        testmin=test
                xmax=0
                for pro in prog:
                    test=numpy.max(pro[:,0])
                    if(test>xmax):
                        xmax=test
                ymax=0
                for pro in prog:
                    test=numpy.max(pro[:,1])
                    if(test>ymax):
                        ymax=test
                if(xmax==0):
                    print('No plotable data')
                    return
            except:
                print('Error in setting axes')
                return

            jet=cm.get_cmap('jet')

            lines=[]
            zed=[]
            for pro in prog:
                pro=numpy.array(pro)
                for i in range(len(pro)-1):
                    lines.append((list(pro[i:i+2,0]),list(pro[i:i+2,1])))
                    zed.append(pro[i+1,2])

            zed=numpy.array(zed)
            lines = [numpy.column_stack([x, y]) for x, y in lines]
            
            #lines=[zip(x,y) for x,y in lines]

            lines=LineCollection(lines,array=zed,cmap=jet)
            self.ax.add_collection(lines)

            #if(self.ref=='subgraphMode' or self.ref=='polishMode'):
            self.ax.set_title("%s" % (title),fontsize=8)

            try:
                self.ax.set_xscale('log')
            except:
                pass
            self.ax.set_xlim(1,xmax*1.05)
            self.ax.set_ylim(-2,ymax+1)

            cb=self.fig.colorbar(lines)
            cb.set_label("NMR EdgeScore",fontsize=8)
            self.ax.set_xlabel("IterationNumber",fontsize=8)
            self.ax.set_ylabel("G1node",fontsize=8)

            self.fig.text(0.15,0.15,"current best score: %i" % testmax,fontsize=8)
            #plt.text(0.02, 0.5, textstr, fontsize=14, transform=plt.gcf().transFigure)

        elif(seleRbox==1):
            in1,in2=self.GetConv()
            if(in1==0):
                print('No optimisation data')
                return
            jet=cm.get_cmap('jet')
            #print(in1,in2)

            norm = matplotlib.colors.Normalize(vmin=0.0, vmax=in1[len(in1)-1][0])#normalise colour bar


            in1=numpy.array(in1)
            in2=numpy.array(in2)
            for i in range(len(in2)):
                ii=len(in2)-i-1
                self.ax.bar(in2[ii][:,1],in2[ii][:,2],color=(jet(norm(1+in2[ii][0,0]))),alpha=1)

            y = numpy.array([1,1+len(in1)])
            colors = cm.jet(y / float(max(y)))
            sm = cm.ScalarMappable(cmap=cm.jet, norm=matplotlib.colors.Normalize(vmin=1, vmax=(1+len(in1))))
            sm._A = []
            self.fig.colorbar(sm)


            self.ax2=self.ax.twinx()
            self.ax2.set_ylabel('Ratio',color='r')


            cax=self.ax2.plot(in1[:,1],in1[:,0],color='r',label='best',linewidth=1)


            if(self.ref=='subgraphMode' or self.ref=='polishMode'):
                self.ax.set_title("optimisation: %s" % (self.subgraph.GetStringSelection()),fontsize=8)
            else:
                self.ax.set_title("optimisation ",fontsize=8)


            self.ax.set_xlabel("Size of MCES",fontsize=8)
            self.ax.set_ylabel("Count",fontsize=8)
            self.ax2.set_ylabel('Largest MCES progress',color='r')

            inbig=in1[len(in1)-1,1]
            for i in range(len(in1)):
                if(in1[i,1]==inbig):
                    imin=i
                    break

            self.fig.text(0.15,0.75,"largest mces found: %i" % inbig,fontsize=8)
            self.fig.text(0.15,0.73, "discovered in iteration: %i" % imin,fontsize=8)

        elif(seleRbox==4):

            self.fig.clear()

            seq,pred,phi,psi,R2,QH,QE,QL,conf,ss=self.ReadTalos()


            pos=range(len(seq))
            self.axA = self.fig.add_subplot(321)
            self.axA.set_axis_off()
            self.axA.text(0.1,7.5,'Prediction Accuracy',fontsize=8)
            self.axA.text(0.1,3.5,'Secondary Structure',fontsize=8)

            for i,se in enumerate(seq):
                rpos=(i-0.5,1)
                if(ss[i]=='L'):
                    col='red'
                elif(ss[i]=='E'):
                    col='blue'
                else:
                    col='green'

                rec1 = plt.Rectangle(rpos,1,2,color=col)
                self.axA.add_artist(rec1)

            for i,se in enumerate(seq):
                rpos=(i-0.5,5)
                print(pred[i])
                if(pred[i]=='Good'):
                    col='green'
                elif(pred[i]=='Warn'):
                    col='red'
                elif(pred[i]=='Dyn'):
                    col='blue'
                else:
                    col='white'

                rec1 = plt.Rectangle(rpos,1,2,color=col)
                self.axA.add_artist(rec1)

            self.axA.set_xlim(-0.5,len(pos)+0.5)
            self.axA.set_ylim(0,10)

            self.axB = self.fig.add_subplot(323)
            self.axB.clear()
            self.axB.plot(pos,QH,color='g')
            self.axB.plot(pos,QE,color='b')
            self.axB.plot(pos,QL,color='r')
            self.axB.set_xlim(-0.5,len(pos)+0.5)
            self.axB.set_ylabel("ss",fontsize=8)


            self.axC = self.fig.add_subplot(325)
            self.axC.plot(pos,R2)
            self.axC.set_xlim(-0.5,len(pos)+0.5)

            self.axC.set_xlabel('Residue',fontsize=8)
            self.axC.set_ylabel('R2',fontsize=8)


            self.ax2 = self.fig.add_subplot(122)
            self.ax2.scatter(phi,psi)
            self.ax2.set_xlabel('phi',fontsize=8)
            self.ax2.set_ylabel('psi',fontsize=8)
            self.ax2.set_xlim(-180,180)
            self.ax2.set_ylim(-180,180)





        self.canvas.draw()

    def ReadTalos(self):
        talos=[]
        pred=[]
        phi=[]
        psi=[]
        QH=[]
        QE=[]
        QL=[]
        R2=[]
        ss=[]
        conf=[]
        seq=[]
        infile='pred.tab'

        if(os.path.exists(infile)==0): #execute talos
            outy=open('talosScript.csh','w')
            outy.write('talos+ -offset -in results.res.peaks.tab')
            outy.close()
            os.system('csh talosScript.csh')

        cnt=0
        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                if(cnt==1):
                    if(len(test)==11):
                        seq.append(test[1])
                        pred.append(test[10])
                        phi.append(float(test[2]))
                        psi.append(float(test[3]))
                        R2.append(float(test[7]))
                if(test[0]=='FORMAT'):
                    cnt+=1
        inny.close()

        cnt=0
        infile='predSS.tab'
        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                if(cnt==1):
                    if(len(test)==9):
                        QH.append(float(test[4]))
                        QE.append(float(test[5]))
                        QL.append(float(test[6]))
                        conf.append(float(test[7]))
                        ss.append(test[8])
                if(test[0]=='FORMAT'):
                    cnt+=1
        phi=numpy.array(phi)
        psi=numpy.array(psi)
        QH=numpy.array(QH)
        QE=numpy.array(QE)
        QL=numpy.array(QL)
        R2=numpy.array(R2)
        return seq,pred,phi,psi,R2,QH,QE,QL,conf,ss

    def ReadProg(self,infile):
        dat=[]
        dit=[]
        inny=open(infile)
        for line in inny.readlines():
            test=line.split()
            if(len(test)==0):
                if(len(dit)!=0):
                    dat.append(numpy.array(dit))
                dit=[]

            if(len(test)>0):
                x=float(test[0])
                y=float(test[2])
                c=float(test[4])
                dit.append((x,y,c))


        if(len(dat)==0):
            return numpy.array(dit)
        else:
            if(len(dit)!=0):
                dat.append(numpy.array(dit))
            return dat

    #get progress data
    def GetProg(self):

        dat=[]
        if(self.ref=='subgraphMode' or self.ref=='polishMode' or self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()
            if(sele=='all'):
                for i in range(len(self.parent.inst.subgraphRef.keys())):
                    ii=len(self.parent.inst.subgraphRef.keys())-i-1

                    #print('looking for subgraph',ii+1,' and nproc',int(self.parent.processors.GetValue()))
                    if(int(self.parent.processors.GetValue())==1):
                        testfile=self.testdir+'/mces_%i.txt.G' % ((ii+1))
                        if(os.path.exists(testfile)):
                            #print('reading ',testfile)
                            dat.append(self.ReadProg(testfile))
                            return dat,testfile
                    else:
                        testfile=self.testdir+'/mces_%i.txt.G.%i' % ((ii+1),0)
                        if(os.path.exists(testfile)): #if we are still running...
                            for j in range((int(self.parent.processors.GetValue()))):
                                #print('reading ',testfile)
                                testfile=self.testdir+'/mces_%i.txt.G.%i' % ((ii+1),j)
                                dat.append(self.ReadProg(testfile))
                            return dat,testfile
                        testfile=self.testdir+'/mces_%i.txt.G' % ((ii+1),)
                        if(os.path.exists(testfile)):
                            #print('reading',testfile)
                            dit=self.ReadProg(testfile)
                            try:
                                #print(dit.shape)
                                dat.append(dit)
                                return dat,testfile
                            except:
                                return dit,testfile
            else:
                #for i in range(len(self.parent.inst.subgraphRef.keys())):
                    #print('looking for subgraph',sele,' and nproc',int(self.parent.processors.GetValue()))
                    if(int(self.parent.processors.GetValue())==1):
                        testfile=self.testdir+'/mces_%s.txt.G' % ((sele))
                        if(os.path.exists(testfile)):
                            #print('reading ',testfile)
                            dat.append(self.ReadProg(testfile))
                            return dat,testfile
                    else:
                        testfile=self.testdir+'/mces_%s.txt.G.%i' % ((sele),0)
                        if(os.path.exists(testfile)): #if we are still running...
                            for j in range((int(self.parent.processors.GetValue()))):
                                #print('reading ',testfile)
                                testfile=self.testdir+'/mces_%s.txt.G.%i' % ((sele),j)
                                dat.append(self.ReadProg(testfile))
                            return dat,testfile
                        testfile=self.testdir+'/mces_%s.txt.G' % ((sele),)
                        if(os.path.exists(testfile)):
                            print('reading',testfile)
                            dit=self.ReadProg(testfile)
                            if(len(dit)==int(self.parent.processors.GetValue())):
                                return dit,testfile
                            else:
                                dat.append(dit)
                                return dat,testfile



        else:
            if(int(self.parent.processors.GetValue())==1):
                testfile=self.testdir+'/mces.txt.G'
                if(os.path.exists(testfile)):
                    dat.append(self.ReadProg(testfile))
                return dat,testfile
            else:
                testfile=self.testdir+'/mces.txt.G.%i' % (0)
                if(os.path.exists(testfile)):
                    for j in range((int(self.parent.processors.GetValue()))):
                        testfile=self.testdir+'/mces.txt.G.%i' % (j)
                        if(os.path.exists(testfile)):
                            dat.append(self.ReadProg(testfile))
                    return dat,testfile
                testfile=self.testdir+'/mces.txt.G'
                if(os.path.exists(testfile)):
                    return self.ReadProg(testfile),testfile

        self.GetLastFile() #if the last file is in community mode...
        if(len(self.fil.split('community'))>1):
            testfile=self.testdir+'/'+self.fil
            return self.ReadProg(testfile),testfile


    def PopulateNOE(self):
        result_dict={}
        self.noe=[]
        if(self.ref=='subgraphMode' or self.ref=='polishMode' or self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()

            if(sele=='all'):
                noe=[]
                for ii in range(len(self.parent.inst.subgraphRef.keys())):
                    testfile=self.testdir+'/G2hist_'+str(ii+1)+'.hist'
                    if(os.path.exists(testfile)):
                        self.GetDistHistNOE(testfile)
            else:
                testfile=self.testdir+'/G2hist_'+str(sele)+'.hist'
                if(os.path.exists(testfile)):
                    self.GetDistHistNOE(testfile)

        else:
            testfile=self.testdir+'/G2hist.hist'
            if(os.path.exists(testfile)):
                self.GetDistHistNOE(testfile)

        #reset listbox
        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)

        cnt=0
        for i,noe in enumerate(self.noe):
            cnt+=1
            num_items = self.source.GetItemCount()
            self.source.InsertStringItem(num_items,str(cnt))
            self.source.SetStringItem(num_items,0,str(noe[0])) #id
            self.source.SetStringItem(num_items,1,str(noe[1])) #A
            self.source.SetStringItem(num_items,2,str(noe[2])) #B
            self.source.SetStringItem(num_items,3,str(noe[7])) #intensity
            self.source.SetStringItem(num_items,4,str(noe[4])) # number
            self.source.SetStringItem(num_items,5,str(noe[6])) #w/s
            self.source.SetStringItem(num_items,6,'%.2f' % (noe[5]*1./(1.*self.nsoln)) ) #hits
            self.source.SetStringItem(num_items,7,str(noe[9])) #dist
            self.source.SetStringItem(num_items,8,str(noe[10])) #stdev

            if(int(noe[5])==0):
                color = (int(255),0, 0)
                self.source.SetItemBackgroundColour(num_items,color)

            """
            self.source.InsertColumn(0, 'id', width = 50,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(1, 'A', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(2, 'B', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(3, 'intensity', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(4, '#', width = 50,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(5, 'tp', width = 50,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(6, 'hits', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(7, 'dist', width = 80,format=wx.LIST_FORMAT_CENTRE)
            self.source.InsertColumn(8, 'stdev', width = 80,format=wx.LIST_FORMAT_CENTRE)
            """
            #stry=''
            #for val in vals:
            #    stry+=val+' '
            #self.source.SetStringItem(num_items,3,stry)

        return

    def PopulateFilt(self):
        result_dict={}
        if(self.ref=='subgraphMode' or self.ref=='polishMode' or self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()

            if(sele=='all'):
                for ii in range(len(self.parent.inst.subgraphRef.keys())):
                    testfile=self.testdir+'/filter_'+str(ii+1)+'.res'
                    if(os.path.exists(testfile)):
                        analysis.AddFile(result_dict,testfile)
            else:
                testfile=self.testdir+'/filter_'+str(sele)+'.res'
                if(os.path.exists(testfile)):
                    analysis.AddFile(result_dict,testfile)
        else:
            testfile=self.testdir+'/filter.Full.res'
            if(os.path.exists(testfile)):
                analysis.AddFile(result_dict,testfile)

        #reset listbox
        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)

        cnt=0
        for key in sorted(result_dict,key=lambda k:len(result_dict[k])):
            vals=result_dict[key]
            cnt+=1
            num_items = self.source.GetItemCount()
            self.source.InsertStringItem(num_items,str(cnt))
            self.source.SetStringItem(num_items,0,str(cnt))
            self.source.SetStringItem(num_items,1,key)
            self.source.SetStringItem(num_items,2,'filter')
            stry=''
            for val in vals:
                stry+=val+' '
            self.source.SetStringItem(num_items,3,stry)

    def UpdateScores(self):
        tig=0  #set to 1 if all or not subgraph/polish modes
        if( self.ref=='subgraphMode' or self.ref=='polishMode' ):
            sele=self.subgraph.GetStringSelection()
            if(sele=='all'):
                tig=1
        elif(self.ref=='nudgeMode'):
            sele=self.subgraph.GetStringSelection()
            if(sele=='all'):
                tig=1
        else:
            tig=1

        ret=0 #return flag
        if(tig==1): #if not subgraph/polish nor sele='all'
            self.result_dict={}
            self.result_dict,nsoln,edges,edgestot,edgesG2=analysis.AnalAll(self.testdir)
        else:
            self.result_dict={}
            mcesfile=self.testdir+'/mces_%s.txt' % (int(sele))
            vf2file=self.testdir+'/vf2_%s.txt' % (int(sele))

            ret=0
            if(os.path.exists(mcesfile)):
                nsoln,edges,edgestot,edgesG2,g1node=analysis.AddFile(self.result_dict,mcesfile)
                mode='mces'
            elif(os.path.exists(vf2file)):
                mode='vf2 '
                nsoln,edges,edgestot,edgesG2,g1node=analysis.AddFile(self.result_dict,vf2file)
            else:
                print('Neither vf2 nor mces file exists for subgraph',sele)
                self.GetLastFile()
                nsoln=0
                edges=0
                edgestot=0
                edgesG2=0
                if(len(self.fil.split('vf2'))>1):
                    mode='vf2'
                else:
                    mode='mces'
                ret=1

        #print(len(self.result_dict.keys()))
        self.conf=0
        self.tot=len(self.result_dict.keys())
        for key,vals in self.result_dict.items():
            if(len(vals)==1):
                self.conf+=1

        self.nsoln=nsoln
        self.confLbl.SetLabel("Confident:     %i / %i" % (self.conf,self.tot))
        self.solnLbl.SetLabel("Solutions:     %i" % nsoln)
        self.nmreLbl.SetLabel("NMR EdgeScore: %i / %i" % (edges,edgestot))
        self.pdbeLbl.SetLabel("PDB EdgeScore: %i" % edgesG2)

        #get current result_dict
        if(len(self.result_dict.keys())==0):
            ret=1

        return ret




    def PopulateList(self,shift='n'):
        #reset listbox
        num_items = self.source.GetItemCount()
        for i in range(num_items):
            self.source.DeleteItem(0)

        if(self.UpdateScores()):
            return

        #populate list
        resultConf={}
        conf=0
        result_dict_temp = copy.deepcopy(self.result_dict)
        for key,vals in self.result_dict.items():
            if(len(vals)==1):
                conf+=1
                resultConf[key]=vals
                del(result_dict_temp[key])
        self.result_dict = result_dict_temp
        #print('conf',conf)
        cnt=0
        for i in range(self.colVal+1):
            if(os.path.exists(self.progress[i][1]+'/confident.res')): #is this mode done?
                #print('testing:',self.progress[i][1]+'/confident.res')
                new={}
                inny=open(self.progress[i][1]+'/confident.res')
                for line in inny.readlines():
                    test=line.split(':')
                    key=test[0].split()[0]
                    ass=test[1].split()[0]
                    if(key in resultConf.keys()):
                        if(ass==resultConf[key][0]):
                            new[key]=ass
                            del(resultConf[key])
                        #else:
                            #print('ODD? confident solution not in final results')
                            #print(key,ass)
                inny.close()

                for key,vals in new.items():
                    cnt+=1

                    num_items = self.source.GetItemCount()
                    if(self.WXV==4):
                        self.source.InsertItem(num_items,str(cnt))
                        self.source.SetItem(num_items,0,str(cnt))
                        self.source.SetItem(num_items,1,key)
                        self.source.SetItem(num_items,2,self.progress[i][0])
                    else:
                        self.source.InsertStringItem(num_items,str(cnt))
                        self.source.SetStringItem(num_items,0,str(cnt))
                        self.source.SetStringItem(num_items,1,key)
                        self.source.SetStringItem(num_items,2,self.progress[i][0])



                    if(shift=='n'):
                        stry=vals
                    else:
                        stry=vals+"(%.2f)" % (self.shiftDict[key][vals]['sc'])
                    if(self.WXV==4):
                        self.source.SetItem(num_items,3,stry)
                    else:
                        self.source.SetStringItem(num_items,3,stry)

                    if(self.progress[i][0]=='subgraphMode' or self.progress[i][0]=='polishMode' or self.progress[i][0]=='nudgeMode' or self.progress[i][0]=='finalMode'):

                        #color = (int(255*percentage), 0, int(255*(1.-percentage)))
                        color = (0,int(255), 0)
                        #self.source.SetCellBackgroundColour(num_items, 2, color)
                        self.source.SetItemBackgroundColour(num_items,color)


                    #    color='green'
                    #else:
                    #    color='yellow'
                    #outlat.write('%s \\begin{footnotesize}\\colorbox{%s!30}{%s}\\end{footnotesize} ' % (vals,color,self.progress[i][0]))
                    #v1=(re.findall(r'\d+',vals)[0])
                    #v2=(re.findall(r'\d+',key)[0])
                    #if(v1!=v2):
                    #    outlat.write('\\begin{footnotesize}\\colorbox{red!30}{%s}\\end{footnotesize} ' % ('INCORRECT!'))
                    #outlat.write('\n\n')

        for key in sorted(self.result_dict,key=lambda k:len(self.result_dict[k])):
            vals=self.result_dict[key]
            cnt+=1
            num_items = self.source.GetItemCount()
            if(self.WXV==4):
                self.source.InsertItem(num_items,str(cnt))
                self.source.SetItem(num_items,0,str(cnt))
                self.source.SetItem(num_items,1,key)
            else:
                self.source.InsertStringItem(num_items,str(cnt))
                self.source.SetStringItem(num_items,0,str(cnt))
                self.source.SetStringItem(num_items,1,key)



            #self.source.SetItem(num_items,2,self.progress[i][0])
            stry=''
            for val in vals:
                if(shift=='n'):
                    stry+=val+' '
                else:
                    stry+=val+"(%.2f)" % (self.shiftDict[key][val]['sc'])
            if(self.WXV==4):
                self.source.SetItem(num_items,3,stry)
            else:
                self.source.SetStringItem(num_items,3,stry)

            #self.source.ForceRefresh()
        #print(len(self.result_dict.keys()))



"""
class NOEMan(wx.App):
    def __init__(self,inherit):
        self.frame_NOEManFrame=NOEManFrame(None,10,'Assignment',inherit)
        self.frame_NOEManFrame.Show(True)
#        return Frame1(parent)

# assign ID numbers
[wxID_FRAME1, wxID_FRAME1BUTTON1, wxID_FRAME1BUTTON2, wxID_FRAME1LISTBOX1,
] = [wx.NewId() for _init_ctrls in range(4)]


class NOEManFrame(wx.Frame):
#    title = 'AssBox'
    def __init__(self,parent, id, title,inherit):
        self._init_ctrls(parent,inherit)
        self.parent=inherit

    def

    def _init_ctrls(self,prnt,parent):
        # BOA generated methods
        wx.Frame.__init__(self, id=wxID_FRAME1, name='', parent=prnt,
              pos=wx.Point(358, 184), size=wx.Size(800, 280),
              style=wx.DEFAULT_FRAME_STYLE, title=u'ListBox Test ...')
        self.SetClientSize(wx.Size(900, 280))


        panel=wx.Panel(self,-1)

        self.ReadG2Hist()

        conn_data=parent.tabOne.conn_data
        corr=[]
        self.corrDict={}
        for i in range(len(conn_data)):
            add=[]

            distppm=conn_data[i].distppm
            conf=1
            if(distppm<0.5):
                conf+=2
            if(distppm<0.1):
                conf+=3

            add.append(conn_data[i].p1)
            add.append(conn_data[i].p2)

            add.append(conn_data[i].s1)
            add.append(conn_data[i].s2)
            add.append(conn_data[i].f1)
            add.append(conn_data[i].f2)
            add.append(conn_data[i].f3)
            add.append(conn_data[i].frac)
            add.append(conn_data[i].distScore)
            add.append(conn_data[i].distppm)
            add.append(conf)
            add.append(0.0)


            self.corrDict[i]=add




        self.lc=SortedListCtrl(panel,self.corrDict)

        self.lc.InsertColumn(0, 'Starting resonance')
        self.lc.InsertColumn(1, 'Ending resonance')
        self.lc.InsertColumn(2, 's/n 1')
        self.lc.InsertColumn(3, 's/n 2')

        self.lc.InsertColumn(4, 'f1(ppm)')
        self.lc.InsertColumn(5, 'f2(ppm)')
        self.lc.InsertColumn(6, 'f3(ppm)')

        self.lc.InsertColumn(7, 'diff')
        self.lc.InsertColumn(8, 'IntDist')
        self.lc.InsertColumn(9, 'Shift')
        self.lc.InsertColumn(10, 'Confidence')
        self.lc.InsertColumn(11, 'Kilter')

        self.lc.SetColumnWidth(0, 140)
        self.lc.SetColumnWidth(1, 153)

        outy=open('confNOE.list','w')

        items=self.corrDict.items()
        for key,data in items:
            num_items = self.lc.GetItemCount()
            self.lc.InsertItem(num_items,str(data[0]))  #add assignment
            self.lc.SetItem(num_items, 0,str(data[0])) #add atom
            self.lc.SetItem(num_items, 1,str(data[1]))  #add atom

            for i in range(len(data)-2):
                self.lc.SetItem(num_items, i+2,str(data[i+2]))  #add atom
            self.lc.SetItemData(num_items, key)



            outy.write('%s\t' % str(data[0]))
            outy.write('%s\t' % str(data[1]))
            for i in range(len(data)-2):
                outy.write('%s\t' % str(data[i+2]))  #add atom
            outy.write('\n')

        outy.close()

        self.Addbutton =  wx.Button(panel, 10, 'Show',(710,10))
        self.Removebutton= wx.Button(panel, 11, 'Remove',(710,60))
        self.Clearbutton = wx.Button(panel, 12, 'Clear',(710,110))
        self.Closebutton = wx.Button(panel, 13, 'Close',(710,160))
        self.Savebutton = wx.Button(panel, 14, 'Save',(710,210))


        #self.pdbfile = infile
        self.textbox = wx.TextCtrl(
            panel,
            size=(150,-1),
            style=wx.TE_PROCESS_ENTER,pos=(690,240))
        #        self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
        self.textbox.SetValue("out/MyList.out")
        #self.textbox.SetValue(self.pdbfile)




        #self.lc.Bind(wx.EVT_LIST_ITEM_SELECTED, self.onItemSelected)
        self.Bind(wx.EVT_LISTBOX_DCLICK, self.OnAdd,self.lc)

        self.Bind (wx.EVT_BUTTON, self.OnAdd, self.Addbutton)
        self.Bind (wx.EVT_BUTTON, self.OnRemove, self.Removebutton)
        self.Bind (wx.EVT_BUTTON, self.OnClear, self.Clearbutton)
        self.Bind (wx.EVT_BUTTON, self.OnClose, self.Closebutton)
        self.Bind (wx.EVT_BUTTON, self.OnSave, self.Savebutton)

        #self.vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.lc, 1, wx.EXPAND)

        vbox=wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.Addbutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Removebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Clearbutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Closebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.Savebutton, 0, wx.ALIGN_CENTER| wx.TOP)
        vbox.Add(self.textbox, 0, wx.ALIGN_CENTER| wx.TOP)
        hbox.Add(vbox)
        panel.SetSizer(hbox)



        #self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        #panel1.SetSizer(self.vbox)

        #hbox  = wx.BoxSizer(wx.HORIZONTAL)
        #hbox.Add(self.Addbutton, 1, wx.EXPAND)
        #hbox.Add(self.Removebutton, 1, wx.EXPAND)
        #hbox.Add(self.Clearbutton, 1, wx.EXPAND)
        #hbox.Add(self.Closebutton, 1, wx.EXPAND)
        #hbox.Add(self.Savebutton, 1, wx.EXPAND)
        #self.vbox.Add(hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        #panel1.SetSizer(self.vbox)

        self.Centre()
        self.Show(True)





        #hbox.Add(vbox1, 1, wx.EXPAND)
        #hbox.Add(vbox2, 1, wx.EXPAND)
        #self.SetSizer(hbox)


#        self.SetSizer(self.vbox)
#        self.vbox.Fit(self)

    def onItemSelected(self, event):
        """"""

        #currentItem = event.m_itemIndex
        #car = self.corrDict[currentItem]
        #print(car)

        #count = self.lc.GetItemCount()
        #self.sorted_artists = [self.list.GetItem(itemId=row, col=0).GetText() for row in xrange(count)]
        #print(self.sorted_artists)
        #print(self.sorted_artists[currentItem])


    def AtoI(self,val):
        for i in range(len(self.parent.tabOne.peak)):
            if(val==self.parent.tabOne.peak[i][0]):
                return i

    def OnAdd(self, event):
        sele=self.lc.GetFirstSelected()
        count = self.lc.GetItemCount()
        col1 = [self.lc.GetItem(itemId=row, col=0).GetText() for row in xrange(count)][sele]
        col2 = [self.lc.GetItem(itemId=row, col=1).GetText() for row in xrange(count)][sele]
        print(col1,col2,self.AtoI(col1))
        self.parent.ComboBox1.SetSelection(self.AtoI(col1))
        self.parent.ComboBox2.SetSelection(self.AtoI(col2))
        #self.parent.NOE=1

        self.parent.on_draw_button(True)


    def OnRemove(self, event):
        index = self.lc.GetFocusedItem()
        self.lc.DeleteItem(index)

    def OnClose(self, event):
        self.Close()

    def OnClear(self, event):
        self.lc.DeleteAllItems()

    def OnSave(self, event):
        self.outfile=self.textbox.GetValue()
        print
        print('Saving list to ',self.outfile)
        outy=open(self.outfile,'w')
        tmpconn=self.parent.tabOne.conn_data
        count=self.lc.GetItemCount()
        for itemIndex in range(count):
#            print(str(self.lc.GetItemText(self.lc.GetItem(itemIndex,0))))
#            print(itemIndex,self.lc.GetItemData(self.lc.GetItem(itemIndex,0)))
            item1 = [self.lc.GetItem(itemId=row, col=0).GetText() for row in xrange(count)][itemIndex]
            item2 = [self.lc.GetItem(itemId=row, col=1).GetText() for row in xrange(count)][itemIndex]
            tick=0
            for i in range(len(tmpconn)):
                if(tmpconn[i].p1==item1 and tmpconn[i].p2==item2):
                    tick=1
                    break
            if(tick==1):
                print('saving entry')
                entry=tmpconn[i]
                outy.write('%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (entry.p1,entry.p2,entry.f1,entry.f2,entry.f3,entry.s1,entry.s2,entry.frac,entry.distScore,entry.distppm,entry.Intscore))
        outy.close()

#            print(itemIndex,self.lc.GetItemText(itemIndex))

#        print(self.lc.GetItem(0))
#        self.lc.DeleteAllItems()



 #           age = item.GetText()





#        self.button2 = wx.Button(id=wxID_FRAME1BUTTON2, label=u'Clear',
#              name='button2', parent=self, pos=wx.Point(104, 312),
#              size=wx.Size(87, 28), style=0)
#        self.button2.Bind(wx.EVT_BUTTON, self.OnButton2Button,
#              id=wxID_FRAME1BUTTON2)

#    def __init__(self, parent):


        return
    def OnListBox1Listbox(self, event):
        '''
        click list item and display the selected string in frame's title
        '''
#        selName = self.listBox1.GetStringSelection()
#        self.SetTitle(selName)
        return

    def OnButton2Button(self, event):
        '''
        click button to clear the listbox items
        '''
        self.listBox1.Clear()
"""
