# from tkinter import CENTER
import wx,string,copy,math,numpy,os
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
import matplotlib.colors as colors
import nmrglue as ng
from matplotlib.figure import Figure
from magma import Magma
# import wx.lib.scrolledpanel as scrolled
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, ArrowStyle
from matplotlib.colors import ListedColormap


class spinSystem(wx.Panel):

    def __init__(self,parent):

        wx.Panel.__init__(self, parent=parent, name="Spin System")

        self.parent=parent.tabMag
        self.inst=Magma(self.parent.magmaParFile,run='n') #get instance of magma
        self.parent=parent.tabOne.molecule
        # self.specs = parent.tabOne.spec
        self.place1=[]
        self.node1 = []
        self.subgraphs = self.inst.subgraphRef.items()
        self.system_graph = self.subgrapher()
        self.create_main_panel()
        self.draw_figure()
        #print(self.subgraphs)
        #sys.exit(100)


    def create_main_panel(self):
        # self.fig=Figure((5, 4), 75)
        # self.canvas = FigCanvas(self, -1, self.fig)

        bs = wx.BoxSizer(wx.VERTICAL)
        self.bs_1 = wx.BoxSizer(wx.HORIZONTAL)

        system_choices = []
        for i in range(len(self.system_graph)):
            system_choices.append(str(i+1)+" ("+str(len(self.system_graph[i]))+")")
        #system_choices = list(numpy.arange(1, len(system_graph)+1))

        self.system_chooser = wx.ComboBox(self, value = system_choices[0], choices = system_choices)
        self.system_chooser.Bind(wx.EVT_COMBOBOX, self.onChoice)

        self.buttonReport = wx.Button(self, label="Report")
        self.buttonReport.Bind(wx.EVT_BUTTON, self.OnButtonReport)

        self.spin_system_chooser_sizer = wx.BoxSizer(wx.VERTICAL)
        self.spin_system_label = wx.StaticText(self, label="Spin system:")
        self.spin_system_chooser_sizer.Add(self.spin_system_label, 1, wx.LEFT | wx.TOP, border=5)
        self.spin_system_chooser_sizer.Add(self.system_chooser,1,  wx.LEFT, border=5)

        # self.assignment_text = wx.StaticText(self, label='Assignment:')
        self.assignment_box = wx.CheckBox(self, label = 'Assignment')
        self.Bind(wx.EVT_CHECKBOX, self.on_assignment_box, self.assignment_box)
        self.spin_system_chooser_sizer.Add(self.assignment_box,1, wx.LEFT| wx.TOP, border=5)

        self.bs_1.Add(self.spin_system_chooser_sizer, 1, wx.LEFT | wx.TOP, border=5)
        self.bs_1.Add(self.buttonReport, 1, wx.LEFT | wx.TOP, border=5)

        self.system_scroll=wx.ScrolledWindow(self,-1)
        self.system_panel=systemPanel(self.system_scroll, self.system_graph)
        width,height = self.system_panel.GetSize()
        unit = 20
        self.system_scroll.SetVirtualSize((width, height))
        self.system_scroll.SetScrollRate(unit, unit)

        # self.system_scroll_horizontal_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # self.system_scroll_horizontal_sizer.Add(self.system_scroll, 1, wx.GROW)
        self.bs_1.Add(self.system_scroll, 1, wx.LEFT | wx.RIGHT | wx.GROW, border=25)
        bs.Add(self.bs_1, 0, wx.LEFT | wx.GROW)


        self.walk_scroll=wx.ScrolledWindow(self,-1, style=wx.VSCROLL)
        # self.walk_scroll = wx.lib.scrolledpanel.ScrolledPanel(self, -1,style=wx.SIMPLE_BORDER)
        # self.walk_scroll.SetupScrolling(scroll_x=False)
        self.walk_panel=walkPanel(self.walk_scroll, self.system_graph, self.parent)

        # self.walk_scroll.EnableScrolling(xScrolling=False, yScrolling=True)
        # self.walk_panel.DoGetBestSize()

        # self.walk_scroll_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # self.walk_scroll_sizer.Add(self.walk_panel, 1, wx.GROW)
        # self.walk_scroll.SetSizer(self.walk_scroll_sizer)


        self.walk_scroll_horizontal_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.walk_scroll_horizontal_sizer.Add(self.walk_scroll, 1, wx.GROW)
        bs.Add(self.walk_scroll_horizontal_sizer, 2, wx.GROW)
        self.SetSizer(bs)

        width,height = self.walk_panel.GetSize()
        print(height)
        # height=200*46
        self.walk_scroll.SetVirtualSize((width+unit, height+unit))
        self.walk_scroll.SetScrollRate(unit, unit)
        bs.Fit(self)
        self.bs_1.Fit(self)
        #self.scroller = walkingPanel(self)
        #self.vbox.Add(self.scroller, 1, wx.LEFT)

        #self.toolbar = NavigationToolbar(self.canvas)


    def OnButtonReport(self,event):
        pass

    def read_in_assignment(self):
        # try:
            self.assignment_file = 'dat/Assigned_interactive.list'
            self.assignment = {}
            for line in open(self.assignment_file, 'r'):
                fields = line.split(' : ')
                peak = fields[0]
                assignment = fields[1].rstrip().split(' ')
                self.assignment[peak] = assignment
            print(self.assignment)
            self.walk_panel.draw_figure(assignment=self.assignment)

        # except:
            # self.assignment_box.SetValue(0)

    def on_assignment_box(self, event):
        if self.assignment_box.IsChecked():
            self.read_in_assignment()
            self.savepdf()
        else:
            self.walk_panel.draw_figure()
            


    def draw_figure(self):
        ""

    def onChoice(self, event):
        a = int(self.system_chooser.GetString(self.system_chooser.GetSelection()).split(' ')[0])-1
        self.system_panel.subgraph_number = a
        self.walk_panel.subgraph_number = a
        #print(self.system_panel.number)


        self.system_panel.draw_figure()
        self.walk_panel.draw_figure()
        width,height = self.walk_panel.GetSize()
        unit = 20
        self.walk_scroll.SetVirtualSize((width+unit, height+unit))
        self.walk_scroll.SetScrollRate(unit, unit)
        self.walk_scroll.Refresh()
        self.walk_scroll.Update()
        self.Refresh()
        self.Update()
        self.savepdf()

    def savepdf(self):
        with PdfPages('foo.pdf') as pdf:
            pdf.savefig(self.walk_panel.fig)



    def GetNodeNumber(self,test,nodeList):
        for j,node in enumerate(nodeList):
            if(node==test):
                return j
        print('could not find node')
        return -1

    def WalkBack(self,sub_noe_node_list,sub_noe_adjacency,Ass=False):
        node=sub_noe_node_list[0] #take this as first guess
        nodeVal=0                 #take this as first guess
        cnt=0
        while(1==1):
            cnt+=1
            tick=0
            print(cnt,node,nodeVal)
            if(cnt==50):
                pass
                #sys.exit(100)
            for adj in sub_noe_adjacency[nodeVal]:
                if(self.inst.NMR.noes[node][adj][2]=='b'):
                    node=adj
                    nodeVal=self.GetNodeNumber(node,sub_noe_node_list)
                    tick=1
                    break
            if(tick==0):
                break
        return node,nodeVal

    def subgrapher(self):
        self.inst.subgraphRef.items()
        cnt = 0
        xMax = 0
        yrun= 0
        results = []
        for key,vals in self.inst.subgraphRef.items(): #for each subgraph...

            sub_noe_node_list=vals['nodes']
            sub_noe_adjacency=vals['adj']
            results_local = []

            node,nodeVal=self.WalkBack(sub_noe_node_list,sub_noe_adjacency)

            print('start node:',node)
            cenMain=numpy.array((0,cnt))

            #wire up the circles.
            jobs={}
            jobs[0]=[]
            jobs[0].append((node,nodeVal,0))
            ylevel=0
            xlevel={}
            xlevel[0]=0
            radMain=1
            done=0
            noddy=0
            while(1==1):
                noddy+=1
                if(noddy==100):
                    break
                node,nodeVal,xpos=jobs[ylevel][0] #unpack current place in the queue

                if(xpos>xMax):
                    xMax=xpos
                jobs[ylevel].pop(0) #remove entry zero from front of queue

                xlevel[ylevel]=xpos

                cenMain=(xpos,ylevel+yrun)
                col='r'
                lob=node

                results_local.append([cenMain,col,lob,radMain]) #draw circle and label
                self.place1.append(cenMain)
                self.node1.append(node)

                cnt=0
                for adj in sub_noe_adjacency[nodeVal]:

                    if(self.inst.NMR.noes[node][adj][2]=='f'):

                        for j in range(len(sub_noe_node_list)):
                            if(adj==sub_noe_node_list[j]):
                                newVal=j
                                break
                        if(cnt not in jobs.keys()):
                            jobs[cnt]=[]

                        if(cnt==ylevel):
                            jobs[cnt].append((adj,newVal,xlevel[ylevel]+2))
                        else:
                            jobs[cnt].append((adj,newVal,xlevel[ylevel]))
                        cnt+=1

                xlevel[ylevel]+=2
                jump=0
                while(1==1):
                    jump+=1

                    if(len(jobs.keys())==0): #if no jobs, we're done.
                        done=1
                        break

                    if(ylevel in list(jobs.keys())): #is current level in jobs?
                        if(len(jobs[ylevel])==0):
                            del jobs[ylevel]

                    if(len(jobs.keys())==0): #if no jobs, we're done.
                        done=1
                        break
                    else:
                        ylevel=numpy.min(list(jobs.keys()))
                    if(jump==100):
                        break

                if(done==1):
                    break
            results.append(results_local)

        yrun+=4
        return results

class systemPanel(wx.Panel):

    def __init__(self,parent, system_graph):

        wx.Panel.__init__(self, parent, -1, size=(15000, 100))
        self.system_graph = system_graph
        self.max_size = 0

        # for i in range(1):
        #     print(len(self.system_graph[0]))
        #     if self.max_size < len(self.system_graph[0]):
        self.max_size = len(self.system_graph[0])
        # size = self.max_size*96.
        self.fig = Figure(figsize=(self.max_size/1.3,96./96.), dpi = 96)
        self.canvas = FigCanvas(self, -1, self.fig)

        self.subgraph_number = 0

        self.vbox = wx.BoxSizer(wx.HORIZONTAL)


        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.sizer = self.vbox

        self.axes = self.fig.add_subplot(111)


        self.draw_figure()



    def AddCircle(self,axes,cen,col,lab,rad):
        circle1 = plt.Circle(cen, rad, color='coral')
        axes.add_artist(circle1)
        axes.text(cen[0],cen[1],lab,fontsize=10,horizontalalignment='center',verticalalignment='center', color='white')


    def draw_figure(self):
        # self.fig.clf()
        # self.max_size = len(self.system_graph[self.subgraph_number])
        # self.fig = Figure(figsize=(self.max_size/96.,96./96.), dpi = 96)

        self.axes.clear()
        self.axes.axis('off')
        # self.axes.set_xlim(-len(self.system_graph[0])/2.-2,len(self.system_graph[0])+2)
        self.axes.set_xlim(-1.5,2*self.max_size-0.5)
        self.axes.set_ylim(-1.2,1.2)
        self.axes.set_facecolor('green')

        #xMax = len(self.system_graph[self.subgraph_number])
        for result in self.system_graph[self.subgraph_number]:

            self.AddCircle(self.axes,result[0],result[1],result[2],result[3]) #draw circle and label


        # self.fig.subplots_adjust(left=(1.0/xMax), right = 1-(1.0/xMax), wspace=0.01, bottom=0.01, top = 0.99)
        #self.size = (100,100)
        self.fig.tight_layout()

        self.canvas.draw()
        self.vbox.Fit(self)



        # self.SetupScrolling()

# class walkScroller(wx.ScrolledWindow):


class walkPanel(wx.Panel):
    def find_max(self, spectrum, full_peak):
        N_coord = numpy.argmin(numpy.abs(self.molecule.spec[spectrum].index1-full_peak[spectrum][0].f2))
        H_coord = numpy.argmin(numpy.abs(self.molecule.spec[spectrum].index2-full_peak[spectrum][0].f1))
        C_coord = numpy.argmax(numpy.fabs(self.molecule.spec[spectrum].data[:, N_coord, H_coord]))
        return self.molecule.spec[spectrum].data[C_coord, N_coord, H_coord]

    def plotter(self, axis,axis_above, spectrum, full_peak, color, factor=False, offset = 0.0):

        if(len(full_peak[spectrum])==0):
            return

        #print(full_peak)
        N_coord = numpy.argmin(numpy.abs(self.molecule.spec[spectrum].index1-full_peak[spectrum][0].f2))
        H_coord = numpy.argmin(numpy.abs(self.molecule.spec[spectrum].index2-full_peak[spectrum][0].f1))
        if factor !=False:

            factor = abs(factor/numpy.max(numpy.fabs(self.molecule.spec[spectrum].data[:, N_coord, H_coord])))
        else:
            factor = 1
            # 
        axis.plot(self.molecule.spec[spectrum].index0-offset, self.molecule.spec[spectrum].data[:, N_coord, H_coord]*factor, color = color, lw=0.5, label=spectrum.upper())

        length_of_chain = len(self.system_graph[self.subgraph_number])
        # print(length_of_chain)

        axis.spines['right'].set_visible(False)
        axis.spines['left'].set_visible(False)
        axis.spines['top'].set_visible(False)
        width = 0.05
        # print('hnca: ', full_peak[spectrum])
        for pk3 in full_peak[spectrum]:
            C_coord = numpy.argmin(numpy.abs(self.molecule.spec[spectrum].index0-pk3.f3))
            inty = self.molecule.spec[spectrum].data[C_coord, N_coord, H_coord]*factor
            different_specs = ['hncacb', 'hncocacb', 'hncanh', 'hncocanh']
            if(spectrum not in different_specs):
                if pk3.tp == 'main':
                    
                    # if pk3.name[-1:] != 'a':  ## Prevent peaks added from other spectra being plotted
                    #     axis.bar(pk3.f3,pk3.inty*factor, width=width, color=color)
                        # axis.scatter(pk3.f3,pk3.inty*factor, marker='x', color=color)

                    if axis_above != 0:
                        xy1 = (pk3.f3,max(0,pk3.inty*factor))
                        ymin, ymax = axis.get_ylim()
                        xym1 = (pk3.f3, axis_above.get_ylim()[0])
                        con = ConnectionPatch(xyA=xym1, xyB=xy1, coordsA="data", coordsB="data",
                            axesA=axis_above, axesB=axis, color=color, lw=1,  facecolor=None, zorder = 100000,  ls=(0,(5,0)), arrowstyle=ArrowStyle("simple", head_length=.6, head_width=.6, tail_width=.05))
                        # axis_above.add_artist(con)
                        axis.add_artist(con)
                        axis.set_ylim(ymin, ymax)

            if(spectrum=='hncacb'):
                if(pk3.tp=='PosMax'):
                    if axis_above != 0:                       
                        xy1 = (pk3.f3,max(0,pk3.inty*factor))
                        ymin, ymax = axis.get_ylim()
                        xym1 = (pk3.f3, axis_above.get_ylim()[0])
                        con = ConnectionPatch(xyA=xym1, xyB=xy1, coordsA="data", coordsB="data",
                        axesA=axis_above, axesB=axis, color=color, lw=1,  facecolor=None, zorder = 100000,  ls=(0,(5,0)), arrowstyle=ArrowStyle("simple", head_length=.6, head_width=.6, tail_width=.05))
                        # axis_above.add_artist(con)
                        axis.add_artist(con)
                        axis.set_ylim(ymin, ymax)
            
            

            if(spectrum=='hncanh'):
                if(pk3.tp=='diag'):
                    if axis_above != 0:     
                        ymin, ymax = axis.get_ylim()                  
                        xy1 = (pk3.f3,ymax)
                        xym1 = (pk3.f3, axis_above.get_ylim()[0])
                        con = ConnectionPatch(xyA=xym1, xyB=xy1, coordsA="data", coordsB="data",
                        axesA=axis_above, axesB=axis, color='black', lw=0.5,  facecolor=None, zorder = 100000,  ls=(0,(5,0)), arrowstyle=ArrowStyle("simple", head_length=.6, head_width=.6, tail_width=.05))
                        # axis_above.add_artist(con)
                        axis.add_artist(con)
                        axis.set_ylim(ymin, ymax)
                
                if(pk3.tp=='plus'):
                    if axis_above != 0:     
                        ymin, ymax = axis.get_ylim()                  
                        xy1 = (pk3.f3,max(0,pk3.inty*factor))
                        xym1 = (pk3.f3, axis_above.get_ylim()[0])
                        con = ConnectionPatch(xyA=xym1, xyB=xy1, coordsA="data", coordsB="data",
                        axesA=axis_above, axesB=axis, color='black', lw=0.5,  facecolor=None, zorder = 100000,  ls=(0,(5,0)), arrowstyle=ArrowStyle("simple", head_length=.6, head_width=.6, tail_width=.05))
                        # axis_above.add_artist(con)
                        axis.add_artist(con)
                        axis.set_ylim(ymin, ymax)
            

                    
            
            # elif(spectrum=='hncacb'):
            #     if(pk3.tp == 'PosMin'):
            #         if axis_above != 0:
            #             xy1 = (pk3.f3,max(0,pk3.inty*factor))
            #             ymin, ymax = axis.get_ylim()
            #             xym1 = (pk3.f3, axis_above.get_ylim()[0])
            #             con = ConnectionPatch(xyA=xym1, xyB=xy1, coordsA="data", coordsB="data",
            #                 axesA=axis_above, axesB=axis, color=color, lw=1,  facecolor=None, zorder = 100000,  ls=(0,(5,0)), arrowstyle=ArrowStyle("simple", head_length=.6, head_width=.6, tail_width=.05))
            #             # axis_above.add_artist(con)
            #             axis.add_artist(con)
            #             axis.set_ylim(ymin, ymax)


            # else:
                # if pk3.name[-1:] != 'a': ## Prevent peaks added from other spectra being plotted
                #     axis.bar(pk3.f3,pk3.inty*factor, width=width, color=color)
                    # axis.scatter(pk3.f3,pk3.inty*factor, marker='x', color=color)

                # axis.text(pk3.f3, (pk3.inty*factor)/4., 'Peak i-1', rotation = 90, va='center', ha='right')
            # print(pk3.name, spectrum, N_coord, pk3.f2, pk3.tp, pk3.f3, pk3.inty, factor, self.find_max(spectrum, full_peak))

    def __init__(self,parent, system_graph, molecule):
        self.parent = parent
        self.system_graph = system_graph
        self.subgraph_number = 0
        # self.spec = spec
        size = 0
        for i in range(len(self.system_graph)):
            # print(len(self.system_graph[i]))
            if size < len(self.system_graph[i]):
                size = len(self.system_graph[i])

        length = 200.*size

        wx.Panel.__init__(self, parent, -1, size=(1200, 400))

        self.fig = Figure(figsize=(1200./96.,length/96.), dpi = 96)
        self.canvas = FigCanvas(self, -1, self.fig)
        self.canvas.mpl_connect('button_release_event', self.on_click)
        # self.toolbar = NavigationToolbar(self.canvas)
        self.vbox = wx.BoxSizer(wx.HORIZONTAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP )
        
        if('hncocacb' and 'hncacb' in molecule.spec.keys()):
            self.cacb = True
        else:
            self.cacb = False
        
        if('hncanh' and 'hncocanh' in molecule.spec.keys()):
            self.HNHN = True
        else:
            self.HNHN = False
        
        self.sizer = self.vbox
        self.SetSizer(self.vbox)
        print('system_graph length', len(self.system_graph[self.subgraph_number]))



        self.molecule = molecule
        self.bore_array = []
        self.draw_figure()


    def on_click(self, event):

        x = event.xdata
        # print(event.y)
        if x != None:
            if event.xdata < 100.0:
                c_coord = numpy.argmin(numpy.abs(self.molecule.spec['hnca'].index0-x))
            else:
                c_coord = numpy.argmin(numpy.abs(self.molecule.spec['hncaco'].index0-x))
            win = TestPopup(self.GetTopLevelParent(), wx.SIMPLE_BORDER, self.molecule, c_coord, x)

        pos = (event.x, event.y)
        pos = (0,0)
        sz =  0
        win.Position(pos, (0, sz))

        win.Show(True)

    # def fold_carbons(self, value, spectrum):
    #     if spectrum == 'hnco':
    #         spectrum = 'hncaco'
        

    def plot_assignment(self, axis, spectrum, assignment, height, cmap):

        def gradient_image(ax, extent, direction=1, cmap_range=(0, 1), **kwargs):
            phi = direction * numpy.pi / 2
            v = numpy.array([numpy.cos(phi), numpy.sin(phi)])
            X = numpy.array([[v @ [1, 0], v @ [1, 1]],
                        [v @ [0, 0], v @ [0, 1]]])
            a, b = cmap_range
            X = a + (b - a) / X.max() * X
            im = ax.imshow(X, extent=extent, interpolation='bicubic',
                        vmin=0, vmax=1, **kwargs)
            # return im

        def gradient_bar(ax, x, y, cmap, width=5, bottom=0, direction = 1):
            for left, top in zip(x, y):
                right = left + width
                print(left, right, top, bottom)
                gradient_image(ax, extent=(left, right, bottom, top),
                        cmap=cmap, cmap_range=(0, 1.0), direction = direction)

        three_letter = self.molecule.p1to3[assignment[0][-1:]]
        bmrb = self.molecule.bmrbC
        carbons = bmrb[three_letter]
        CA = carbons['CA']
        if 'CB' in carbons.keys():
           CB = carbons['CB']
        CO = carbons['C']


        
        my_cmap = cmap(numpy.arange(cmap.N))
        my_cmap[:,-1] = 0.2 #numpy.linspace(0, 1, cmap.N)
        my_cmap = numpy.vstack((my_cmap, my_cmap[::-1]))
        my_cmap_for = ListedColormap(my_cmap)
        # my_cmap_bak = ListedColormap(my_cmap[::-1])
        
        # print(CO, CA, CB)
        if spectrum[-2:] == 'ca':
            folded = self.molecule.DoUnAlias_number(CA[0], 'hnca')
            # axis.axvline(folded-(CA[1]), lw= 1, ls='--', color='k')
            # axis.axvline(folded+(CA[1]), lw= 1, ls='--', color='k')
            gradient_bar(axis, [folded-(CA[1]),], [height,], width = CA[1]*2., cmap=my_cmap_for, direction = 1)
            
            axis.set_aspect('auto')
        if spectrum[-2:] == 'co':
            
            folded=self.molecule.DoUnAlias_number(CO[0], 'hncaco')
            # axis.axvline(folded-(CO[1]), lw= 1, ls='--', color='k')
            # axis.axvline(folded+(CO[1]), lw= 1, ls='--', color='k')
            gradient_bar(axis, [folded-(CO[1]),], [height,], width = CO[1]*2., cmap=my_cmap_for, direction = 1)
            
            axis.set_aspect('auto')
            

    def draw_figure(self, assignment=False):

        self.fig.clf()
        size = len(self.system_graph[self.subgraph_number])*200.
        self.SetSize((1200,size))
        self.fig.set_size_inches(1200./96., size/96., forward=True)
        if not self.cacb or 'hncacb' not in self.molecule.spec.keys() or 'hncocacb' not in self.molecule.spec.keys():
            self.axes = self.fig.subplots(len(self.system_graph[self.subgraph_number]), 2, sharex='col')
        elif(self.HNHN==False):
            self.axes = self.fig.subplots(len(self.system_graph[self.subgraph_number]), 3, sharex='col')
            self.axes[0,2].set_title('Alpha and Beta Carbons')
        else:
            self.axes = self.fig.subplots(len(self.system_graph[self.subgraph_number]), 4, sharex='col')
            self.axes[0,2].set_title('Alpha and Beta Carbons')
            self.axes[0,3].set_title('Backbone Nitrogens')


        self.axes[0,0].set_title('Alpha Carbons', fontsize=14)
        self.axes[0,1].set_title('Carbonyl Carbons', fontsize=14)
        self.popup = wx.PopupWindow(self)
        for j, peak in enumerate(list(self.system_graph[self.subgraph_number])):
            name = peak[2]
            j = len(self.system_graph[self.subgraph_number])-j-1
            for i, pk in enumerate(self.molecule.peak.keys()):
                if peak[2] == pk:
                    full_peak = self.molecule.peak[pk]
                    
                    try:
                        if assignment == False:
                            self.axes[j,0].set_ylabel(full_peak['hnco'][0].name[:-2], fontsize=14)
                        elif len(assignment[name]) == 1:
                            self.axes[j,0].set_ylabel(full_peak['hnco'][0].name[:-2]+' ('+assignment[name][0]+')', fontsize = 14)
                        self.axes[j,0].axis('on')
                        self.axes[j,0].set_yticks([])

                        self.axes[j,1].set_yticks([])
                    except:
                        import sys
                        print('Ensure that any manual peaks that have been added have been inputted into all 2D and 3D peaklists ')
                        #sys.exit(100)

                    if self.cacb:
                        self.axes[j,2].set_yticks([])
                        self.axes[j,2].axis('on')
                    hnca_max_val = 0.0
                    hnco_max_val = 0.0
                    hncacb_max_val = 0.0
                    hncanh_max_val = 0.0

                    boring_array = []
                    if 'hnca' in full_peak:
                        if j >0:
                            self.plotter(self.axes[j,0], self.axes[j-1,0], 'hnca', full_peak, 'darkblue')
                            
                            
                            
                        else:
                            self.plotter(self.axes[j,0], 0, 'hnca', full_peak, 'darkblue')
                        hnca_max_val = self.find_max('hnca', full_peak)
                        # if assignment!=False:
                        #         print(name)
                        #         self.plot_assignment(self.axes[j,0], 'hnca', assignment[name], hnca_max_val, cmap=cm.Blues)
                        
                    if 'hncoca' in full_peak:
                        
                        if j >0:
                            if numpy.fabs(hnca_max_val) > 0.:
                                self.plotter(self.axes[j,0], self.axes[j-1,0], 'hncoca', full_peak, 'red', hnca_max_val, self.molecule.HNCAmed)
                            else:
                                self.plotter(self.axes[j,0], self.axes[j-1,0], 'hncoca', full_peak, 'red')
                        else:
                            if numpy.fabs(hnca_max_val) > 0.:
                                self.plotter(self.axes[j,0], 0, 'hncoca', full_peak, 'red', hnca_max_val, self.molecule.HNCAmed)
                            else:
                                self.plotter(self.axes[j,0], 0, 'hncoca', full_peak, 'red')
                        


                    if 'hnco' in full_peak:
                        if j > 0:
                            self.plotter(self.axes[j,1], self.axes[j-1,1], 'hnco', full_peak, 'g')
                            
                            
                            
                        else:
                            
                            self.plotter(self.axes[j,1], 0, 'hnco', full_peak, 'g')
                        hnco_max_val = self.find_max('hnco', full_peak)
                        # if assignment!=False:
                        #             self.plot_assignment(self.axes[j,1], 'hnco', assignment[name], hnco_max_val, cmap=cm.Purples)
                        


                    if 'hncaco' in full_peak:
                        if j > 0:
                            print()
                            if numpy.fabs(hnco_max_val) > 0.:
                                self.plotter(self.axes[j,1], self.axes[j-1,1], 'hncaco', full_peak, 'purple', hnco_max_val, self.molecule.HNCOmed)
                                
                            else:
                                self.plotter(self.axes[j,1], self.axes[j-1,1], 'hncaco', full_peak, 'purple')
                        else:
                            if numpy.fabs(hnco_max_val) > 0.:
                                self.plotter(self.axes[j,1], 0, 'hncaco', full_peak, 'purple', hnco_max_val, self.molecule.HNCOmed)
                            else:
                                self.plotter(self.axes[j,1], 0, 'hncaco', full_peak, 'purple')
                        


                    if self.cacb and 'hncacb' in full_peak:
                        if j>0:
                            self.plotter(self.axes[j,2], self.axes[j-1,2], 'hncacb', full_peak, 'darkblue')
                            hncacb_max_val = self.find_max('hncacb', full_peak)
                        else:
                            self.plotter(self.axes[j,2], 0, 'hncacb', full_peak, 'darkblue')
                        

                        

                    if self.cacb and 'hncocacb' in full_peak:
                        if j>0:
                            if numpy.fabs(hncacb_max_val) > 0.:
                                self.plotter(self.axes[j,2], self.axes[j-1,2], 'hncocacb', full_peak, 'darkorange', hncacb_max_val, self.molecule.HNCACBmed)
                            else:
                                self.plotter(self.axes[j,2], self.axes[j-1,2], 'hncocacb', full_peak, 'darkorange')
                        else:
                            self.plotter(self.axes[j,2], 0, 'hncocacb', full_peak, 'darkorange')
                        



                    
                    if self.HNHN and 'hncanh' in full_peak:
                        self.axes[j,3].set_yticks([])
                        if j > 0:
                            try:
                                self.axes[j,3].axvline(self.molecule.peak[pk]['hncanh'][0].f2, color='orchid',linewidth=3.0)

                            except:
                                pass
                            self.plotter(self.axes[j,3], self.axes[j-1,3], 'hncanh', full_peak, 'darkblue')
                            hncanh_max_val = self.find_max('hncanh', full_peak)
                            
                        else:
                            try:
                                self.axes[j,3].axvline(self.molecule.peak[pk]['hncanh'][0].f2, color='orchid', linewidth=3.0)

                            except:
                                pass
                            self.plotter(self.axes[j,3], 0, 'hncanh', full_peak, 'darkblue')
                            



                    if self.HNHN and 'hncocanh' in full_peak:
                        if j > 0:
                            if numpy.fabs(hncanh_max_val) > 0.:
                                self.plotter(self.axes[j,3], self.axes[j-1,3], 'hncocanh', full_peak, 'darkgreen', hncanh_max_val, self.molecule.HNCANHmed)
                            else:
                                self.plotter(self.axes[j,3], self.axes[j-1,3], 'hncocanh', full_peak, 'darkgreen')
                        else:
                            self.plotter(self.axes[j,3], 0, 'hncocanh', full_peak, 'darkgreen')
                        
 

        xhigh, xlow = self.axes[0,0].get_xlim()[::-1]
        self.axes[0,0].set_xlim(xhigh, xlow)
        self.axes[0,0].legend(frameon=False)
        self.axes[0,1].legend(frameon=False)
        if(self.cacb==True):
            self.axes[0,2].legend(frameon=False)
        if(self.cacb==True and self.HNHN==True):
            self.axes[0,3].legend(frameon=False)

        xhigh, xlow = self.axes[0,1].get_xlim()[::-1]
        self.axes[0,1].set_xlim(xhigh, xlow)
        self.axes[0,0].autoscale(enable=True, axis='x', tight=True)
        # self.axes[0,0].autoscale(enable=True, axis='y', tight=True)
        self.axes[0,1].autoscale(enable=True, axis='x', tight=True)
        # self.axes[0,1].autoscale(enable=True, axis='y', tight=True)
        if(self.cacb==True and self.HNHN==True):
            self.axes[0,2].autoscale(enable=True, axis='x', tight=True)
            self.axes[0,3].autoscale(enable=True, axis='x', tight=True)

        self.axes[len(self.system_graph[self.subgraph_number])-1,0].set_xlabel(r'$\delta_{C\alpha}$/ppm')
        self.axes[len(self.system_graph[self.subgraph_number])-1,1].set_xlabel(r'$\delta_{CO}$/ppm')
        if(self.cacb):
            self.axes[len(self.system_graph[self.subgraph_number])-1,2].set_xlabel(r'$\delta_{C\alpha/C\beta}$/ppm')
        if(self.HNHN):
            self.axes[len(self.system_graph[self.subgraph_number])-1,3].set_xlabel(r'$\delta_{N}$/ppm')
        
        self.fig.tight_layout(h_pad=0.079)
        self.fig.subplots_adjust(hspace=0.08)

        self.canvas.draw()
        # self.vbox.Fit(self)
            #exit()


class TestPopup(wx.PopupWindow):


    def get_text_positions(self, text, x_data, y_data, txt_width, txt_height):
        a = zip(y_data, x_data)
        text_positions = list(y_data)
        for index, (y, x) in enumerate(a):
            local_text_positions = [i for i in a if i[0] > (y - txt_height)
                                and (abs(i[1] - x) < txt_width * 2) and i != (y,x)]
            if local_text_positions:
                sorted_ltp = sorted(local_text_positions)
                if abs(sorted_ltp[0][0] - y) < txt_height: #True == collision
                    differ = numpy.diff(sorted_ltp, axis=0)
                    a[index] = (sorted_ltp[-1][0] + txt_height, a[index][1])
                    text_positions[index] = sorted_ltp[-1][0] + txt_height*1.01
                    for k, (j, m) in enumerate(differ):
                        #j is the vertical distance between words
                        if j > txt_height * 2: #if True then room to fit a word in
                            a[index] = (sorted_ltp[k][0] + txt_height, a[index][1])
                            text_positions[index] = sorted_ltp[k][0] + txt_height
                            break
        return text_positions

    def text_plotter(self, text, x_data, y_data, text_positions, txt_width,txt_height):
        for z,x,y,t in zip(text, x_data, y_data, text_positions):
            self.axes.annotate(str(z), xy=(x-txt_width/2, t), size=12)
            if y != t:
                self.axes.arrow(x, t,0,y-t, color='red',alpha=0.3, width=txt_width*0.1,
                    head_width=txt_width, head_length=txt_height*0.5,
                    zorder=0,length_includes_head=True)


    def GetLevels(self,min_level,fac,ctr_level):
        levels=[]
        levels.append(min_level)
        for i in range(ctr_level-1):
            levels.append(levels[i]*fac)
        levels=numpy.array(levels)
        levels=numpy.concatenate((-1*levels[::-1],levels)) #reflect on negative axis
        return levels
    #----------------------------------------------------------------------
    def __init__(self, parent, style, molecule, c_coord, c_ppm):
        """Constructor"""
        wx.PopupWindow.__init__(self, parent, style)
        self.molecule = molecule
        panel = wx.Panel(self)
        self.panel = panel
        panel.SetBackgroundColour("CADET BLUE")

        self.fig = Figure(figsize=(600./96., 450./96.), dpi=96)
        self.canvas = FigCanvas(self, -1, self.fig)
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.axes = self.fig.subplots(1)
        if c_ppm < 100.0:
            spec = 'hnca'
            noise =  self.molecule.spec[spec].noise
            color = 'r'
        else:
            spec = 'hncaco'
            noise = self.molecule.spec[spec].noise
            color='g'
        Xs=self.molecule.spec[spec].XX_proj
        Ys=self.molecule.spec[spec].YY_proj
        Zs=self.molecule.spec[spec].data[c_coord,:,:]
        x_res =self.molecule.spec[spec].index2[0]-self.molecule.spec[spec].index2[1]
        y_res =self.molecule.spec[spec].index1[0]-self.molecule.spec[spec].index1[1]

        levels = self.GetLevels(noise, 1.2, 3)

        self.axes.contour(Ys, Xs, Zs, levels=levels, colors=color)
        #print(self.molecule.spec['hnca'].peak)
        plotting_text = []
        plotting_x = []
        plotting_y = []
        for peak in self.molecule.spec[spec].peak:
                if numpy.abs(peak.f3p-c_ppm) < 0.3:
                    plotting_text.append(peak.name)
                    plotting_x.append(peak.f1+x_res*2)
                    plotting_y.append(peak.f2+y_res*2)
        text_positions = self.get_text_positions(plotting_text, plotting_x, plotting_y, x_res*10., y_res*10.)
                    # self.axes.text(peak.f1+x_res*2, peak.f2+y_res*2, peak.name)
        self.text_plotter(plotting_text, plotting_x, plotting_y, text_positions, x_res*10., y_res*10.)
        #exit()
        #self.axes.plot(numpy.sin(range(0,100)))
        # sz = st.GetBestSize()
        # self.SetSize( (sz.width+20, sz.height+20) )
        # panel.SetSize( (sz.width+20, sz.height+20) )
        self.SetSize( (600, 450) )
        panel.SetSize( (600, 450) )

        panel.Bind(wx.EVT_LEFT_DOWN, self.OnMouseLeftDown)
        panel.Bind(wx.EVT_MOTION, self.OnMouseMotion)
        panel.Bind(wx.EVT_LEFT_UP, self.OnMouseLeftUp)
        panel.Bind(wx.EVT_RIGHT_UP, self.OnRightUp)

        self.cid1 = self.canvas.Bind(wx.EVT_LEFT_DOWN, self.OnMouseLeftDown)
        self.cid2 = self.canvas.Bind(wx.EVT_MOTION, self.OnMouseMotion)
        self.cid3 = self.canvas.Bind(wx.EVT_LEFT_UP, self.OnMouseLeftUp)
        self.cid4 = self.canvas.Bind(wx.EVT_RIGHT_UP, self.OnRightUp)

        wx.CallAfter(self.Refresh)

    def OnMouseLeftDown(self, evt):
        self.Refresh()
        self.ldPos = evt.GetEventObject().ClientToScreen(evt.GetPosition())
        self.wPos = self.ClientToScreen((0,0))
        self.panel.CaptureMouse()

    def OnMouseMotion(self, evt):
        if evt.Dragging() and evt.LeftIsDown():
            dPos = evt.GetEventObject().ClientToScreen(evt.GetPosition())
            nPos = (self.wPos.x + (dPos.x - self.ldPos.x),
                    self.wPos.y + (dPos.y - self.ldPos.y))
            self.Move(nPos)

    def OnMouseLeftUp(self, evt):
        if self.panel.HasCapture():
            self.panel.ReleaseMouse()

    def OnRightUp(self, evt):

        if self.HasCapture():
            self.ReleaseMouse()
        if self.panel.HasCapture():
            self.panel.ReleaseMouse()

        self.canvas.ReleaseMouse()
        self.Show(False)
        self.Destroy()
