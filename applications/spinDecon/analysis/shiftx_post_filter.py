import linecache
#import pdb
import collections,numpy
import sys
#from mces import GenericMethods
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
#dictNMR = collections.defaultdict(list)
class shiftXNMR():
    def __init__(self, shiftXFile, inputFile, inputChains, inputResidues, combinedResults):
        self.chains = inputChains
        self.residues = inputResidues #residues
        self.shiftXFile = shiftXFile  #shiftX shifts
        self.NMRFile = inputFile      #exp shifts
        self.shiftXDict = {}
        self.NMRDict = {}
        self.CLimit = 1.44 #0.9754
        self.resultsFile = combinedResults  #current assignment file


    def addToDict(self, resTypeShort, resTypeMed, atomName, dicty):
        dicty[resTypeShort] = {}
        dicty[resTypeShort]['resType'] = resTypeMed
        dicty[resTypeShort]['atomName'] = atomName

    def ReadShiftxFile(self):
        #CREATE SHIFT DICTIONARY
        chainID = ""
        for line in open(self.shiftXFile, 'r'):
            if len(line.split(",")) == 4:
                (resNumber, resType, atomName, nmrValue) = line.split(",")
                areChains = False
            elif len(line.split(",")) == 5:
                (chainID, resNumber, resType, atomName, nmrValue) = line.split(",")
                areChains = True
            else:
                continue

            #print chainID,resNumber,resType,atomName,nmrValue
            #if it's one of the residues we care about
            if (resType in list(self.resTypeDict.keys())) and (self.resTypeDict[resType]['resType'] in self.residues):
                if (not chainID):
                    self.populateShiftXDict(self.shiftXDict, self.resTypeDict, resType, atomName, resNumber, nmrValue)
                elif chainID in self.chains:
                    if chainID not in list(self.shiftXDict.keys()):
                        self.shiftXDict[chainID] = {}
                    self.populateShiftXDict(self.shiftXDict[chainID], self.resTypeDict, resType, atomName, resNumber, nmrValue)




    def SetMethylDict(self):
        self.resTypeDict = {} #indected by single letter aa
        self.addToDict('I', 'ILE', ('CD1','HD1'), self.resTypeDict)
        self.addToDict('L', 'LEU', ('CD1', 'CD2','HD1','HD2'), self.resTypeDict)
        self.addToDict('V', 'VAL', ('CG1', 'CG2','HG1','HG2'), self.resTypeDict)
        self.addToDict('A', 'ALA', ('CB','HB'), self.resTypeDict)
        #self.addToDict('T', 'THR', ('CG2',), resTypeDict)
        #self.addToDict('M', 'MET', ('CE',), resTypeDict)
        #resNumbers = []

    def Parse(self):

        self.SetMethylDict()
        #resNumbers = []

        self.ReadShiftxFile()
        """
        #CREATE SHIFT DICTIONARY
        chainID = ""
        for line in open(self.shiftXFile, 'r'):
            if len(line.split(",")) == 4:
                (resNumber, resType, atomName, nmrValue) = line.split(",")
                areChains = False
            elif len(line.split(",")) == 5:
                (chainID, resNumber, resType, atomName, nmrValue) = line.split(",")
                areChains = True
            else:
                continue

            #print chainID,resNumber,resType,atomName,nmrValue
            #if it's one of the residues we care about
            if (resType in resTypeDict.keys()) and (resTypeDict[resType]['resType'] in self.residues):
                if (not chainID):
                    self.populateShiftXDict(self.shiftXDict, resTypeDict, resType, atomName, resNumber, nmrValue)
                elif chainID in self.chains:
                    if chainID not in self.shiftXDict.keys():
                        self.shiftXDict[chainID] = {}
                    self.populateShiftXDict(self.shiftXDict[chainID], resTypeDict, resType, atomName, resNumber, nmrValue)
        """
        #sorted by residuenumber
        for key,vals in list(self.shiftXDict.items()):
            print(key,vals)

        #CREATE NMR DATA DICITONARY
        for line in open(self.NMRFile, 'r'):
            if len(line.split()) == 3:
                (resID, atomName, nmrValue) = line.split()
            else:
                continue
            if resID not in list(self.NMRDict.keys()):
                self.NMRDict[resID] = {}
                self.NMRDict[resID]['resType'] = resID[-1]
                self.NMRDict[resID]['name'] = {}
            self.NMRDict[resID]['name'][atomName]= float(nmrValue)
            #print resID,atomName,nmrValue

        #sorted by residue number
        #for key,vals in self.NMRDict.items():
        #    print key,vals

        #CURRENT ASSIGNMENTS
        resDict = self.readInCombRes()


        #CORRELATE THE TWO
        full={}
        for key,vals in list(resDict.items()):
            if(key not in list(full.keys())):
                full[key]={}
            #print
            #print key,vals
            for val in vals:
                full[key][val]={}
            #print self.NMRDict[key]
            self.correlate(key,vals,'C',full)
            self.correlate(key,vals,'H',full)

        print()
        for key,vals in list(full.items()):
            for val in vals:
                print(key,val)
                if(len(full[key][val]['H'])==2):
                    tp=list(full[key][val]['C'].keys())[0][1] #position label

                    t2='1'  #get score for 1
                    namH1='H'+tp+t2
                    namC1='C'+tp+t2
                    sc1=( (full[key][val]['H'][namH1][2]*10.)**2.+full[key][val]['C'][namC1][2]**2.)**0.5

                    t2='2'  #get score for 2
                    namH2='H'+tp+t2
                    namC2='C'+tp+t2
                    sc2=( (full[key][val]['H'][namH2][2]*10.)**2.+full[key][val]['C'][namC2][2]**2.)**0.5

                    if(sc2>sc1):
                        #print 'sc1 wins'
                        del full[key][val]['C'][namC2]
                        del full[key][val]['H'][namH2]
                        full[key][val]['sc']=sc1
                    else:
                        #print 'sc2 wins'
                        del full[key][val]['C'][namC1]
                        del full[key][val]['H'][namH1]
                        full[key][val]['sc']=sc2

                else:
                    print(list(full[key][val]['C'].keys()))
                    tp=list(full[key][val]['C'].keys())[0][1] #position label

                    if(tp=='B'): #for alanines
                        t2=''
                    else:
                        t2='1'  #get score for 1
                    namH1='H'+tp+t2
                    namC1='C'+tp+t2
                    print(list(full[key][val]['H'].keys()))
                    print(list(full[key][val]['C'].keys()))

                    v=0
                    try:
                        v+=(full[key][val]['H'][namH1][2]*10.)**2.
                    except:
                        pass

                    try:
                        v+=full[key][val]['C'][namC1][2]**2.
                    except:
                        pass

                    sc1=(v)**0.5
                    if(sc1==0):
                        sc1='-'
                    full[key][val]['sc']=sc1


        #print
        #for key,vals in full.items():
            #print key,':'
            #for val in vals:
                #print val,'-',full[key][val]['sc']
                #print full[key][val]['H'][full[key][val]['H'].keys()[0]]
                #print full[key][val]['C'][full[key][val]['C'].keys()[0]]

        return full

    #correlate dictionarys  and calculate differences
    def correlate(self,key,vals,targ,full):
        print('shit')
        for NMRnam in list(self.NMRDict[key]['name'].keys()):
            if(NMRnam[0]==targ):
                sh=self.NMRDict[key]['name'][NMRnam]
                for val in vals:
                    if(targ not in list(full[key][val].keys())):
                        full[key][val][targ]={}

                    #print val
                    print(val)
                    resN=val.split('C')[0]
                    #print self.shiftXDict[resN]
                    print(resN)
                    #if(resN not in self.shiftXDict.keys()):
                    for PDBnam in list(self.shiftXDict[resN]['name'].keys()):
                        if(PDBnam[0]==targ):
                            ph=self.shiftXDict[resN]['name'][PDBnam]
                            #print PDBnam,ph,NMRnam,sh,numpy.fabs(ph-sh)
                            full[key][val][targ][PDBnam]=sh,ph,numpy.fabs(ph-sh)


        return 0


    def readInCombRes(self):
        newResDict = {}
        for line in open(self.resultsFile, 'r'):
            key, vals = line.split(":")
            vals = vals.split()
            keyStrip = key.strip()
            if keyStrip not in list(newResDict.keys()):
                newResDict[keyStrip] = []
            newVals = []
            for val in vals:
                newResDict[keyStrip].append(val)
        return newResDict


    def populateShiftXDict(self, dictS, resTypeDict, resType, atomName, resNumber, nmrValue):
        #if(len(atomName)!=3): #give up if name is too short
        #    return

        resID=str(resNumber)

        if atomName not in resTypeDict[resType]['atomName']:
            return

        if resID not in list(dictS.keys()):
            dictS[resID] = {}
            dictS[resID]['type'] = resType #resTypeDict[resType]['resType'][0]
            dictS[resID]['name'] = {}
        if atomName not in list(dictS[resID].keys()):
            dictS[resID]['name'][atomName] = float(nmrValue)
