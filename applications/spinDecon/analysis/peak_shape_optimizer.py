import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.optimize as opt
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
class PeakShapeOptimizer():
    def __init__(self, spec, XX, YY, ZZ,indexes, peak_locations, starting_values):
        self.spec = spec
        self.values = starting_values

        print(self.values)
        self.XX = XX
        self.YY = YY
        self.ZZ = ZZ
        self.indexes = indexes
        self.unpack_init()
        self.index_s = [np.ceil(self.s[0] / np.abs(self.indexes[0][1] - self.indexes[0][0])),
                        np.ceil(self.s[1] / np.abs(self.indexes[1][1] - self.indexes[1][0])),
                        np.ceil(self.s[2] / np.abs(self.indexes[2][1] - self.indexes[2][0]))]

        i = 0
        for peak_location in peak_locations:

            self.peak_location = peak_location

            self.initialise_peak_space()
            # print self.peak.shape, self.index_s, self.peak_location
            self.values2 = opt.leastsq(self.GetChi, x0=self.values, maxfev=100)[0]
            print(self.values2)
            if i == 0:
                self.x0 = self.values2
            else:
                self.x0 = np.vstack((self.x0, self.values2))
            i+=1
            self.values = self.values2

        self.x0 = np.average(self.x0, axis=0)
        self.x0[0:3] = np.abs(np.sin(self.x0[0:3]))

    def initialise_peak_space(self):
        self.minx = max(0,int(self.peak_location[0] - 5 * self.index_s[0]))
        self.miny = max(0,int(self.peak_location[1] - 5 * self.index_s[1]))
        self.minz = max(0, int(self.peak_location[2] - 5 * self.index_s[2]))
        self.maxx = min(self.spec.shape[0], int(self.peak_location[0] + 5 * self.index_s[0]))
        self.maxy = min(self.spec.shape[1], int(self.peak_location[1] + 5 * self.index_s[1]))
        self.maxz = min(self.spec.shape[2], int(self.peak_location[2] + 5 * self.index_s[2]))

        self.peak_XX = self.XX[self.minx:self.maxx, self.miny:self.maxy, self.minz:self.maxz]
        self.peak_YY = self.YY[self.minx:self.maxx, self.miny:self.maxy, self.minz:self.maxz]
        self.peak_ZZ = self.ZZ[self.minx:self.maxx, self.miny:self.maxy, self.minz:self.maxz]

        self.peak_meshes = np.array([self.peak_XX, self.peak_YY, self.peak_ZZ])

        self.peak = self.spec[self.minx:self.maxx, self.miny:self.maxy, self.minz:self.maxz]

        self.peak_data = np.array([self.peak, self.peak, self.peak])

        self.peak_indexes = np.array([self.indexes[0][self.minx:self.maxx], self.indexes[1][self.miny:self.maxy],
                                      self.indexes[2][self.minz:self.maxz]])

        self.x = np.ones_like(self.peak_XX) * self.indexes[0][self.peak_location[0]]
        self.y = np.ones_like(self.peak_XX) * self.indexes[1][self.peak_location[1]]
        self.z = np.ones_like(self.peak_XX) * self.indexes[2][self.peak_location[2]]
        self.locations = np.array([self.x, self.y, self.z])
        self.unpack(self.values)
        self.count = 0
    def unpack_init(self):
        self.n = self.values[0:3]
        self.s = self.values[3:6]
        self.r = self.values[6:9]
        self.i = self.values[9]

    def unpack(self, x):
        self.n = np.abs(np.sin(x[0:3]))
        self.s = np.array(x[3:6]).astype(float)/2.355
        self.r = x[6:9]
        self.n_array = np.abs(np.sin(np.array([np.ones_like(self.peak_XX)*x[0], np.ones_like(self.peak_YY)*x[1], np.ones_like(self.peak_ZZ)*x[2]])))
        self.s_array = np.array([np.ones_like(self.peak_XX)*x[3], np.ones_like(self.peak_YY)*x[4], np.ones_like(self.peak_ZZ)*x[5]])/2.355
        self.r_array = np.array([np.ones_like(self.peak_XX)*x[6], np.ones_like(self.peak_YY)*x[7], np.ones_like(self.peak_ZZ)*x[8]])
        self.i = np.array([x[9]])


    def Ycalc(self):
        sim_peak = self.i*np.multiply(np.subtract(1.,self.n_array),np.exp(-np.subtract(self.peak_meshes, self.locations) * np.subtract(self.peak_meshes , self.locations)
                                                                   / np.multiply(2., np.power(self.s_array, 2)))) + self.i*np.multiply(self.n_array, np.divide(np.power(np.divide(self.r_array,2.),2.),(np.power(np.subtract(self.peak_meshes,self.locations),2.)+np.power(np.divide(self.r_array,2.),2.))))
        self.ycalc = np.multiply(sim_peak[0,:,:,:], np.multiply(sim_peak[1,:,:,:], sim_peak[2,:,:,:]))




    def GetChi(self, x):
        self.unpack(x)


        if self.count % 100 == 0:
            print(self.count, self.r)
        self.Ycalc()
        # plt.plot(np.sum(np.sum(((self.peak / np.max(self.peak)) - self.ycalc), axis=1), axis=0))
        # plt.plot(np.sum(np.sum((self.ycalc), axis=0), axis=0), ls='--')
        #
        # plt.plot(np.sum(np.sum((self.peak/ np.max(self.peak)), axis=0), axis=0), color="k")
        # plt.show()
        self.count +=1
        return (((self.peak/np.max(self.peak))-self.ycalc)).flatten()
