import wx
import os
"""
Academic Use Licence

These licence terms apply to all licences granted by THE CHANCELLOR, MASTERS AND SCHOLARS OF THE UNIVERSITY OF OXFORD whose administrative offices are at University Offices, Wellington Square, Oxford OX1 2JD, United Kingdom (the "University") for use of UniDecNMR ("the Software") downloaded from the following website: https://github.com/charliebuchanan/UniDecNMR ("the Website")
By downloading the Software through the Source, you (the "Licensee") are confirming that you agree that your use of the Software is subject to these licence terms.

PLEASE READ THESE LICENCE TERMS CAREFULLY BEFORE DOWNLOADING THE SOFTWARE THROUGH THIS WEBSITE.  IF YOU DO NOT AGREE TO THESE LICENCE TERMS YOU SHOULD NOT DOWNLOAD THE SOFTWARE.

THE SOFTWARE IS INTENDED FOR USE BY ACADEMICS CARRYING OUT RESEARCH AND NOT FOR USE BY CONSUMERS OR COMMERCIAL BUSINESSES.

1.	Academic Use Licence
1.1	The Licensee is granted a limited non-exclusive and non-transferable royalty free licence to download and use the Software provided that the Licensee will:
(a)	limit their use of the Software to their own internal academic non-commercial research which is undertaken for the purposes of education or other scholarly use; 
(b)	not use the Software for or on behalf of any third party or to provide a service or integrate all or part of the Software into a product for sale or license to third parties;
(c)	use the Software in accordance with the prevailing instructions and guidance for use given on the Website and comply with procedures on the Website for user identification, authentication and access;
(d)	comply with all applicable laws and regulations with respect to their use of the Software; and 
(e)	ensure that the Copyright Notice "Copyright (c) 2022, University of Oxford" appears prominently wherever the Software is reproduced and on any documents or other material created using the Software.
1.2	The Licensee may only reproduce, modify, transmit or transfer the Software where:
(a)	such reproduction, modification, transmission or transfer is for academic, research or other scholarly use;
(b)	the conditions of this Licence are imposed upon the receiver of the Software or any modified Software;
(c)	all original and modified Source Code is included in any transmitted software program; and
(d)	the Licensee grants the University an irrevocable, indefinite, royalty free, non-exclusive unlimited licence to use and sub-licence any modified Source Code as part of the Software.

1.3	The University reserves the right at any time and without liability or prior notice to the Licensee to revise, modify and replace the functionality and performance of the access to and operation of the Software.
1.4	The Licensee acknowledges and agrees that the University owns all intellectual property rights in the Software.  The Licensee shall not have any right, title or interest in the Software.
1.5	This Licence will terminate immediately and the Licensee will no longer have any right to use the Software or exercise any of the rights granted to the Licensee upon any breach of the conditions in Section 1 of this Licence.

2.	Indemnity and Liability 
2.1	The Licensee shall defend, indemnify and hold harmless the University against any claims, actions, proceedings, losses, damages, expenses and costs (including without limitation court costs and reasonable legal fees) arising out of or in connection with the Licensee's possession or use of the Software, or any breach of these terms by the Licensee. 
2.2	The Software is provided on an 'as is' basis and the Licensee uses the Software at their own risk. No representations, conditions, warranties or other terms of any kind are given in respect of the the Software and all statutory warranties and conditions are excluded to the fullest extent permitted by law. Without affecting the generality of the previous sentences, the University gives no implied or express warranty and makes no representation that the Software or any part of the Software: (a) will enable specific results to be obtained; or (b) meets a particular specification or is comprehensive within its field or that it is error free or will operate without interruption; or (c) is suitable for any particular, or the Licensee's specific purposes. 
2.3	Except in relation to fraud, death or personal injury, the University's liability to the Licensee for any use of the Software, in negligence or arising in any other way out of the subject matter of these licence terms, will not extend to any incidental or consequential damages or losses, or any loss of profits, loss of revenue, loss of data, loss of contracts or opportunity, whether direct or indirect.
2.4	The Licensee hereby irrevocably undertakes to the University not to make any claim against any employee, student, researcher or other individual engaged by the University, being a claim which seeks to enforce against any of them any liability whatsoever in connection with these licence terms or their subject-matter. 

3.	General 
3.1	Severability - If any provision (or part of a provision) of these licence terms is found by any court or administrative body of competent jurisdiction to be invalid, unenforceable or illegal, the other provisions shall remain in force.
3.2	Entire Agreement - These licence terms constitute the whole agreement between the parties and supersede any previous arrangement, understanding or agreement between them relating to the Software. 
3.3	Law and Jurisdiction - These licence terms and any disputes or claims arising out of or in connection with them shall be governed by, and construed in accordance with, the law of England. The Licensee irrevocably submits to the exclusive jurisdiction of the English courts for any dispute or claim that arises out of or in connection with these licence terms.

If you are interested in using the Software commercially, please contact Oxford University Innovation Limited to negotiate a licence. Contact details are enquiries@innovation.ox.ac.uk 

"""
#FGA
class LoadFilePopUp(wx.Frame):
    """A class for the pop-up window that lets you load a results directory into the GUI.
    """
    def __init__(self, parent):
        # wx.Frame.__init__(self, parent, flags=wx.BORDER_SUNKEN)
        wx.Frame.__init__(self, None, wx.ID_ANY, "Load results",
                          style = wx.DEFAULT_FRAME_STYLE)
        self.Centre(direction=wx.BOTH)

        self.listBox = wx.ListCtrl(self, id=wx.ID_ANY, style=wx.LC_REPORT | wx.EXPAND)
        self.listBox.InsertColumn(0, 'Directory', width=wx.LIST_AUTOSIZE)
        self.listBox.InsertColumn(1, 'Size', width=wx.LIST_AUTOSIZE)

        loadButton = wx.Button(self, label="Load", size=wx.DefaultSize)
        deleteButton = wx.Button(self, label="Delete",size=wx.DefaultSize)
        cancelButton = wx.Button(self, label="Cancel",size=wx.DefaultSize)
        loadButton.Bind(wx.EVT_BUTTON, self.LoadResults)
        deleteButton.Bind(wx.EVT_BUTTON,self.deleteResults)
        cancelButton.Bind(wx.EVT_BUTTON,self.cancel)

        mainSizer = wx.BoxSizer(wx.HORIZONTAL)
        buttonSizer = wx.BoxSizer(wx.VERTICAL)
        buttonSizer.Add(loadButton)
        buttonSizer.Add(deleteButton)
        buttonSizer.Add(cancelButton)
        mainSizer.Add(self.listBox, 0, wx.ALL|wx.EXPAND, 5)
        mainSizer.Add(buttonSizer,0,wx.ALL,5)
        self.SetSizerAndFit(mainSizer)

        self.PopulateListCtrl()

    def PopulateListCtrl(self):
        files = os.listdir('.')
        listOfFiles = []
        for f in files:
            if f.startswith("results"):
                listOfFiles.append(f)
        currentResultsSize = self.folderSize('results')
        #print currentResultsSize
        count = 0
        self.index = 0
        for l in listOfFiles:
            lSplit = l.split('results')
            resSize = self.folderSize(l)
            if len(lSplit) > 1:
                if lSplit[1]:
                    if currentResultsSize == resSize:
                        #print l
                        count +=1
                    name = lSplit[1]
                    line = "Line %s" % self.index
                    self.listBox.InsertStringItem(self.index, line)
                    self.listBox.SetStringItem(self.index, 0, name)
                    self.listBox.SetStringItem(self.index, 1, str(resSize))
                    self.index += 1
        if count == 0:
            msg = """ The current results directory has not been saved.
            Are you sure you want to continue?
            """
            dlg = wx.MessageDialog(self, msg, "Warning", wx.YES_NO)
            ans = dlg.ShowModal()
            if ans == wx.ID_YES:
                dlg.Destroy()
                self.Show(True)
            elif ans == wx.ID_NO:
                dlg.Destroy()
                self.Destroy()
        else:
            self.Show(True)

    def LoadResults(self,e):
        itemIdx = self.listBox.GetFirstSelected()
        if itemIdx != -1:
            resDir = self.listBox.GetItemText(itemIdx, col=0)
            path = 'results' + resDir
            if os.path.exists(path) == 1:
                for item in os.listdir('results'):
                    os.system('rm -r results/'+item)
                for item in os.listdir(path):
                    os.system('cp -r ' + path + '/' + item + ' results')
            else:
                print('Sorry, that directory does not exist.')
            self.Destroy()

    def deleteResults(self,e):
        itemIdx = self.listBox.GetFirstSelected()
        if itemIdx != -1:
            resDir = self.listBox.GetItemText(itemIdx, col=0)
            path = 'results' + resDir
            if os.path.exists(path) == 1:
                os.system('rm -r '+ path)
            else:
                print('Sorry, that directory does not exist.')
            self.Destroy()

    def cancel(self,e):
        self.Destroy()

    def folderSize(self, path):
        total = 0
        for entry in os.listdir(path):
            ePath = path+'/'+entry
            if os.path.isfile(ePath):
                total += os.path.getsize(ePath)
            elif os.path.isdir(ePath):
                total += self.folderSize(ePath)
        return total
