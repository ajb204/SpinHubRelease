# INDIANA Back-End
# (C) Gogulan Karunanithy 2019
# University of Oxford

import numpy
import copy, os, sys
import matplotlib.pyplot as plt
import nmrglue as ng
import scipy.constants
from random import randint
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit, fit_report
from scipy.optimize import brentq
from scipy.special import jv, jvp, jn_zeros, jnp_zeros
import random
import time
from matplotlib import cm
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"]})


class readData():
    def __init__(self, data, gradients, big_delta, little_delta):

        # self.dic, self.data = ng.pipe.read(infile) # read in data
        # self.data = numpy.absolute(self.data)
        # self.size = self.data.shape
        # self.uc0 = ng.pipe.make_uc(self.dic, self.data, dim=0)
        # self.uc1 = ng.pipe.make_uc(self.dic, self.data, dim=1) # make universal constants for easy referencing
        # print(self.uc0.ppm_scale(), self.uc1.ppm_scale())

        # self.ppm = []
        # for x in range(int(self.size[1])):
        #     self.ppm.append(float(self.uc1.ppm(x))) # get ppm value for each point in direct dimension

        # self.ppm = numpy.array(self.ppm)

        # if input_script == False:
        #     from vpar import GetPar
        #     self.gzlvl1 = self.retNumpArray(GetPar('raw/',('gzlvl1',)))
        #     self.delta = float(GetPar('raw/', ('gt1',))[0])
        #     self.BigT = self.retNumpArray(GetPar('raw/', ('BigT',)))
        # else:
        #     self.delta, self.BigT,self.gzlvl1 = self.readInput(input_script)

        self.gzlvl1 = gradients
        self.BigT = big_delta
        self.delta = little_delta


        self.pathExists('out')
        self.pathExists('figs')

        self.params = {}
        self.params['delta'] = self.delta
        self.params['gzlvl1'] = self.gzlvl1
        self.params['BigT'] = self.BigT
        # self.params['dac'] = self.dac


    def readInput(self, infile):
        print('getting diffusion  input values from script...')
        with open(infile, 'r') as inny:
            for line in inny:
                line = line.split()
                if len(line)>0 and line[0]!='#':
                    if line[0]=='delta':
                        delta = float(line[1])
                        print('delta = ', delta)
                    elif line[0]=='BigT':
                        BigT = []
                        for i in range(1,len(line)):
                            BigT.append(float(line[i]))
                        BigT = numpy.array(BigT)
                        print('BigT = ', BigT)
                    elif line[0]=='grad':
                        grad = []
                        for i in range(1,len(line)):
                            grad.append(float(line[i]))
                        grad = numpy.array(grad)
                        print('grad = ', grad)
        try:
            return delta,BigT,grad
        except:
            print('it looks as though the input script is not setup correctly')
            sys.exit(100)

    def pathExists(self, folder): #create 'out' and 'figs' folders'
        if (os.path.exists(folder)==0):
            os.system('mkdir ' +folder)

    def retNumpArray(self, array): #return numpy arrays from vpar
        temp = []
        for i in range(len(array)):
            temp.append(float(array[i]))
        return numpy.array(temp)


    def sortData(self): # sort data by gzlvl1 lowest to highest in case order was randomised
        self.sortedData = numpy.zeros([len(self.BigT), len(self.gzlvl1), len(self.ppm)])

        BigT_index = numpy.argsort(self.BigT)
        gzlvl_index = numpy.argsort(self.gzlvl1)

        for i in range(len(self.BigT)):
            for j in range(len(self.gzlvl1)):
                #self.sortedData[i,j,:] = self.data[BigT_index[i]*len(self.gzlvl1)+ gzlvl_index[j],:] #map original data into sorted array
                self.sortedData[i,j,:] = self.data[BigT_index[i] + gzlvl_index[j]*len(self.BigT),:] #map original data into sorted array

        self.BigT = numpy.sort(self.BigT) #put into ascending order
        self.gzlvl1 = numpy.sort(self.gzlvl1)

        self.params['gzlvl1'] = copy.deepcopy(self.gzlvl1)
        self.params['BigT'] = copy.deepcopy(self.BigT)

    # def sortData_Bruker(self):
    #     self.sortedData = numpy.zeros([len(self.BigT), len(self.gzlvl1), len(self.ppm)])

    #     # BigT_index = numpy.argsort(self.BigT)
    #     gzlvl_index = numpy.argsort(self.gzlvl1)

    #     for i in range(len(self.BigT)):
    #         for j in range(len(self.gzlvl1)):
    #             #self.sortedData[i,j,:] = self.data[BigT_index[i]*len(self.gzlvl1)+ gzlvl_index[j],:] #map original data into sorted array
    #             self.sortedData[i,j,:] = self.data[BigT_index[i] + gzlvl_index[j]*len(self.BigT),:] #map original data into sorted array

    #     self.BigT = numpy.sort(self.BigT) #put into ascending order
    #     self.gzlvl1 = numpy.sort(self.gzlvl1)

    #     self.params['gzlvl1'] = copy.deepcopy(self.gzlvl1)
    #     self.params['BigT'] = copy.deepcopy(self.BigT)


    def getNoise(self, minVal, maxVal):
        topNoise = self.uc1(str(maxVal)+'ppm')
        botNoise = self.uc1(str(minVal)+'ppm')
        self.noise = numpy.std(self.sortedData[0, 0, topNoise:botNoise])
        return self.noise


    def integrateSpec(self, minVal, maxVal): #integrate region and place into numpy.array sorted by [BigT, gzlvl1]
        intTop = self.uc1(str(maxVal)+'ppm')
        intBot = self.uc1(str(minVal)+'ppm')
        self.numPoints = intBot - intTop
        self.intDat = numpy.zeros([len(self.BigT), len(self.gzlvl1)])

        for i in range(len(self.BigT)):
            for j in range(len(self.gzlvl1)):
                self.intDat[i,j] = numpy.sum(self.sortedData[i,j,intTop:intBot]) #integrate data over required range #integrate data over required range


###### exchange/restriction equations

class FitDiff():
    def __init__(self, data, params, meth='cell'):
        self.gyH = scipy.constants.value('proton gyromag. ratio')*1e-4
        self.gzlvl1 = params['gzlvl1']
        self.delta = params['delta']
        self.BigT = params['BigT']
        # self.dac = params['dac']
        self.dat = data

        # print(self.dat)
        # exit()

        self.g = self.gzlvl1
        self.q = self.g*self.gyH*self.delta
        self.q2 = self.q**2.

        if meth=='mono':
            if len(self.BigT)>1:
                self.param_List = ['D1', 'r1']
            else:
                self.param_List = ['D1']

        elif meth=='biexp':
            if len(self.BigT)>1:
                self.param_List = ['D1','D2','r1_1','r1_1']
            else:
                self.param_List = ['D1','D2']

        elif(meth=='3pool'):
        
            self.param_List = ['D1','D2','D3','r1_1','r1_2','r1_3','kab','kba','kbc','kcb','rad1','rad2']
        else:
            self.param_List = ['D1','D2','r1_1','r1_2','kex','rad']

        self.mode=meth

    def Bess(self, x, rad, i):
        from scipy.special import jv, jvp, jn_zeros, jnp_zeros
        return rad*x*jvp(1.5, rad*x, n=1) - 0.5*jv(1.5, rad*x) #actual functions we want to get zeros for

    def findZero(self, rad, num, tol):
        from scipy.special import jv, jvp, jn_zeros, jnp_zeros
        from scipy.optimize import brentq
        zeros = jnp_zeros(2, num) # find zeros for our approx function
        zeros = zeros/rad # scale our zeros
        ans = numpy.zeros(num)

        ans = numpy.array([brentq(self.Bess,zeros[i]-tol,zeros[i]+tol, args = (rad, i), maxiter=10000) for i in range(num)]) #use brentq method to find zeros
        return ans


    def chiFunc_rest(self,params, dat, BigT, q2, BigT_list, boot=False):
        vals = self.exchange_func_rest(params, BigT, q2, BigT_list, boot)
        if(self.mode=='mono'):
            return ((vals) - (dat)).flatten()
        else:
            return (numpy.log(vals) - numpy.log(dat)).flatten()


    def calc_diff_vals(self, rad, D, num, delta, BigT):
        intra = copy.deepcopy(BigT)
        roots = self.findZero(rad, num, 1/rad)
        tay_limit = 1e-3
        atten_vals = 0.
        taylor = 0.
        if (D*BigT[-1]/rad**2.0)<tay_limit:
            intra, val = numpy.meshgrid(BigT, roots)
            # val = copy.deepcopy(roots)# roots[i]
            k = val**6.0*D**3.0*delta**2.0*(intra-delta/3) - 0.5*val**8.0*D**4.0*delta**2.0*intra**2.0 + (1./6.)*val**10.0*D**5.0*intra*delta**2.0*(intra**2.0 + 0.5*delta**2.0)  - (1./24.)*val**12.0*D**6.0*(intra**4.0*delta**2.0 + intra**2.0*delta**4.0)
            k += (1./5040.)*val**14.0*D**7.0*(42*delta**2.0*intra**5.0 + 70*intra**3.0*delta**4.0 + 14*intra*delta**6.0 - 2*delta**7.0)
            k -= (1./10080.)*val**16.0*D**8.0*(14*intra**6.0*delta**2.0 + 35*intra**4.0*delta**4.0 + 14*intra**2.0*delta**6.0)
            k += (1./181440.)*val**18.0*D**9.0*(36*intra**7.0*delta**2.0 + 126*intra**5.0*delta**4.0 + 83*intra**3.0*delta**6.0 + 9*intra*delta**8.0 - delta**9.0)
            k += (1./181440.)*val**18.0*D**9.0*(36*intra**7.0*delta**2.0 + 126*intra**5.0*delta**4.0 + 83*intra**3.0*delta**6.0 + 9*intra*delta**8.0 - delta**9.0)
            k -= (1./362880.)*val**20.0*D**10.0*(3*intra**8.0*delta**2.0 + 14*intra**6.0*delta**4.0 + 14*intra**4.0*delta**6.0 + 3*intra**2.0*delta**8.0)
            mroots = copy.deepcopy(val)

        else:
            mBigT, mroots = numpy.meshgrid(BigT, roots)
            lbit = mroots**2.0*D
            k = 2*lbit*delta - 2 + 2*numpy.exp(-(lbit*self.delta))+ 2*numpy.exp(-(lbit*mBigT)) - numpy.exp(-(lbit*(mBigT+self.delta))) - numpy.exp(-(lbit*(mBigT-self.delta)))

        G = mroots**6.0*(rad**2.0*mroots**2.0 - 2)

        final = k/G
        atten_vals = numpy.sum(final, axis=0)

        little_d = (2.0*atten_vals)/(D**2.0*(BigT-self.delta/3.)*self.delta**2.)

        return little_d


    def exchange_func_rest(self, params, BigT, q2, BigT_list, boot=False):
        if(self.mode=='mono'):

            D_out = params['D1']
            r1_in = params['r1']
            pre_exp = params['pre_exp']
            return pre_exp*numpy.exp(-(numpy.array(BigT*(q2*D_out+r1_in)  ))) # basically just monoexpoential decay here


        elif(self.mode=='3pool'):

            # we send BigT and q in as grids this allows us full flexibility
            # BigT_list is a list of all BigT values

            D_out = params['D1']
            D_in_real = params['D2']
            p_out = params['pa']
            kex = params['kex']
            r1_out = params['r1_1']
            r1_in = params['r1_2']
            pre_exp = params['pre_exp']
            rad = params['rad']

            if p_out >0.99999999:
                return pre_exp*(numpy.array(-q2*BigT*D_out)) # basically just monoexpoential decay here

            p_in = 1-p_out

            R_plus = r1_in + r1_out
            R_minus = r1_in - r1_out

            if rad > 0.:
                # this step is quite slow so we will just calculate for all BigT values and then assign
                # as appropriate instead of calculating multiple times to create a grid
                D_temp = numpy.array(self.calc_diff_vals(rad, D_in_real, self.num,self.delta, BigT_list))
                # D_temp is a list of apparent D_in in order of BigT
                if boot:
                    D_in = numpy.zeros((q2.shape[0],q2.shape[1]))
                    for i,vals in enumerate(BigT_list):
                        D_in[BigT == vals] = D_temp[i] # map values into D_grid
                else:
                    D_in = numpy.tile(D_temp, (q2.shape[1],1)).transpose() # simple tiling here
            else:
                D_in = numpy.ones((q2.shape[0],q2.shape[1]))*D_in_real # create grid of same D_point


            #this section is taking eigenvalues and vectors of the 2x2 rate matrix:
            #R= [ -Da q^2  -Ra-kab        +kba      ]
            #   [    +kab            -Db q^2 -Rb-kba]
            #so it can calculate rho(t)=exp(Rt)*rho(0)
            #we want to do something similar for three pool.
            #'sig' is summing over contributions from all pools.
            #for 3 pool we need to do the same for:
            #R= [ -Da-Ra-kab        +kba            0     ]
            #   [    +kab      -Db-Rb-kba-kbc     +kcb    ]
            #   [     0             +kbc      -Dc-Rc-kcb  ]

            #mission 1: verify that the analytical sig looks like the numerical eigenvalue method
            #mission 2: do the same trick with the 3x3 matrix, and decide on sets of fitting parameters
            
            R=numpy.zeros((2,2))
            R[0,0]=-D_in*q2-R_in - p_in*kex
            R[0,1]= p_in*kex
            R[1,0]= p_out*kex
            R[1,1]= -D_out*q2-R_out-p_out*kex

            s0=numpy.array((pa),(pb))
            sigVec=numpy.dot(expm(R*BigT),s0)
            sigN=numpy.sum(sigVec)


            D_plus = D_in + D_out
            D_minus = D_in - D_out

            xi_plus = q2*D_plus + R_plus + kex*(p_out+p_in)
            xi_minus = q2*D_minus + R_minus + kex*(p_out-p_in)

            psi = numpy.sqrt(xi_minus**2.0 + 4*kex**2.0*p_in*p_out)

            lambda_plus = 0.5*(xi_plus + psi)
            lambda_minus = 0.5*(xi_plus - psi)

            sig = pre_exp*(numpy.exp(-lambda_minus*BigT)+(psi - kex + (p_in-p_out)*(q2*D_minus+R_minus))*((numpy.exp(-lambda_plus*BigT)-numpy.exp(-lambda_minus*BigT))/(2.0*psi)))

            print(sig,sigN)

            return sig




        elif(self.mode=='bi'):
            #self.param_List = ['D1','D2','pa','r1_1','r1_2','pre_exp']

            D_out = params['D1']
            D_in  = params['D2']
            Pa  = params['pa']
            r1_in = params['r1_1']
            r1_out = params['r1_2']
            pre_exp = params['pre_exp']
            return pre_exp*(Pa*numpy.exp(-(numpy.array(BigT*(q2*D_in+r1_in)  )))   +(1-Pa) *numpy.exp(-(numpy.array(BigT*(q2*D_out+r1_out)  )))  ) # basically just monoexpoential decay here



        else:
            # we send BigT and q in as grids this allows us full flexibility
            # BigT_list is a list of all BigT values

            D_out = params['D1']
            D_in_real = params['D2']
            p_out = params['pa']
            kex = params['kex']
            r1_out = params['r1_1']
            r1_in = params['r1_2']
            pre_exp = params['pre_exp']
            rad = params['rad']

            if p_out >0.99999999:
                return pre_exp*(numpy.array(-q2*BigT*D_out)) # basically just monoexpoential decay here

            p_in = 1-p_out

            R_plus = r1_in + r1_out
            R_minus = r1_in - r1_out

            if rad > 0.:
                # this step is quite slow so we will just calculate for all BigT values and then assign
                # as appropriate instead of calculating multiple times to create a grid
                D_temp = numpy.array(self.calc_diff_vals(rad, D_in_real, self.num,self.delta, BigT_list))
                # D_temp is a list of apparent D_in in order of BigT
                if boot:
                    D_in = numpy.zeros((q2.shape[0],q2.shape[1]))
                    for i,vals in enumerate(BigT_list):
                        D_in[BigT == vals] = D_temp[i] # map values into D_grid
                else:
                    D_in = numpy.tile(D_temp, (q2.shape[1],1)).transpose() # simple tiling here
            else:
                D_in = numpy.ones((q2.shape[0],q2.shape[1]))*D_in_real # create grid of same D_point

            D_plus = D_in + D_out
            D_minus = D_in - D_out

            xi_plus = q2*D_plus + R_plus + kex*(p_out+p_in)
            xi_minus = q2*D_minus + R_minus + kex*(p_out-p_in)

            psi = numpy.sqrt(xi_minus**2.0 + 4*kex**2.0*p_in*p_out)

            lambda_plus = 0.5*(xi_plus + psi)
            lambda_minus = 0.5*(xi_plus - psi)

            sig = pre_exp*(numpy.exp(-lambda_minus*BigT)+(psi - kex + (p_in-p_out)*(q2*D_minus+R_minus))*((numpy.exp(-lambda_plus*BigT)-numpy.exp(-lambda_minus*BigT))/(2.0*psi)))

            return sig


    def fitFunc_rest(self, boot_num = 0, para_flg = 'n', ncpus=1):
        from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit, fit_report
        from scipy.optimize import leastsq

        params = Parameters()
        meth = 'leastsq'

        if(self.mode=='mono'):
            self.param_List = ['D1','r1','pre_exp']
            params.add('D1', value = 1.86e-5, min = 0.0, vary = True)
            params.add('r1', value = 1.0, min = 0.0, vary = True)
            params.add('pre_exp', value = self.dat[0,0], min = 0.0, vary = True)

        elif(self.mode=='bi'):

            self.param_List = ['D1','D2','pa','r1_1','r1_2','pre_exp']
            params.add('D1', value = 1.86e-5, min = 0.0, vary = True)
            params.add('D2', value = 5.01e-6, min = 0.0, vary = True)
            params.add('pa', value = 0.1, min = 0.0, max = 1.0, vary = True)
            params.add('r1_1', value = 1.0, min = 0.0, vary = True)
            params.add('r1_2', value = 1.0, min = 0.0, vary = True)
            params.add('pre_exp', value = self.dat[0,0], min = 0.0, vary = True)

        else:
            self.param_List = ['D1','D2','pa','kex','r1_1','r1_2','pre_exp','rad']
            self.num = 20
            params.add('D1', value = 1.86e-5, min = 0.0, vary = True)
            params.add('D2', value = 5.01e-6, min = 0.0, vary = True)
            params.add('pa', value = 0.9, min = 0.0, max = 1.0, vary = True)
            params.add('kex', value = 2.0, min = 1e-10, vary = True)
            params.add('r1_1', value = 1.0, min = 0.0, vary = True)
            params.add('r1_2', value = 1.0, min = 0.0, vary = True)
            params.add('pre_exp', value = self.dat[0,0], min = 0.0, vary = True)
            params.add('rad', value = 10.0e-4, min = 0.0,  vary = True)

        self.BigT_grid = numpy.tile(self.BigT, (self.dat.shape[1],1)).transpose() # make a grid of BigT values
        self.q_grid  = numpy.tile(self.q , (self.dat.shape[0], 1)) # make grid of q values
        self.q2_grid = numpy.tile(self.q2, (self.dat.shape[0], 1)) # make grid of q values

        # sys.exit(100)
        # print(self.dat, self.BigT_grid, self.q2_grid, self.BigT)
        # print(params)
        # print(self.dat, self.BigT_grid, self.q2_grid, self.BigT)
        minner = Minimizer(self.chiFunc_rest, params,(self.dat, self.BigT_grid, self.q2_grid, self.BigT,))
        self.result = minner.minimize(method=meth)

        self.write_results('out/results_covarianceErr.txt', self.result)
        # self.plot_results(result, self.BigT, self.dat, plot = 'figs/results_covarianceErr.pdf')


        if boot_num>0:

            boot_params_all = []
            boot_chi_all = []
            if para_flg == 'y':
                import pp
                ppservers=()
                job_server = pp.Server(ncpus = ncpus, ppservers = ppservers, restart=True)
                init = copy.deepcopy(result.params)
                jobs = [(i, job_server.submit(self.boot_run,(result,i,),(),("numpy","lmfit" ))) for i in range(boot_num)]
                for i,job in jobs:
                    boot_params_all.append(job().params)

                job_server.print_stats()
                time.sleep(2)
                job_server.destroy()

            else:

                for i in range(boot_num):
                    boot_all = self.boot_run(result,i)
                    boot_params_all.append(boot_all.params)
                    boot_chi_all.append(boot_all.chisqr)

            vals_av, vals_std, vals_all = self.unpack_boot(boot_params_all)

            if(self.mode=='mono'):
                param_List = ['D1','r1','pre_exp',]
                for item in param_List:
                    print(item, ' average = ', vals_av[item], ' +/- ', vals_std[item])

            else:
                param_List = ['D1','D2','pa','kex','r1_1','r1_2','pre_exp','rad','perm','cell_dens']
                for item in param_List:
                    print(item, ' average = ', vals_av[item], ' +/- ', vals_std[item])

            self.pathExists('bootstrap_results')
            self.write_results_boot('bootstrap_results/boot_results_rad_all.txt', vals_av, vals_std, param_List)

            self.plot_results_boot(vals_av, vals_std ,self.BigT, self.dat, plot='bootstrap_results/boot_final.pdf')

            best_fit = {}
            bootErr={}

            for params in self.param_List:
                inty,av, sigma, fwhm = self.hist_results(vals_all[params],params)
                best_fit[params] = result.params[params].value
                bootErr[params]=sigma

                #take the smallest error out of standard deviation and the sigma of the normal fit.
                #(std might be too wide if there are outliers)
                if(sigma<vals_std[params]):
                    vals_std[params]=sigma


            self.plot_results_boot(result.params, vals_std ,self.BigT, self.dat, plot='figs/final_results.pdf')
            self.write_results_boot('out/results_bootstrapErr.txt', best_fit, vals_std, self.param_List)

            #self.plot_results_boot(result.params, bootErr ,self.BigT, self.dat, plot='figs/final_results.pdf')
            #self.write_results_boot('out/results_bootstrapErr.txt', best_fit, bootErr, self.param_List)



    def hist_results(self,params,param_name):
        binNo = 50
        bMin = numpy.min(params)*0.95
        bMax = numpy.max(params)*1.05
        bBins = numpy.linspace(bMin, bMax, binNo)

        exhist, exedges = numpy.histogram(params,bBins)

        #setting first and last position in the histogram to be zero...
        exhist[0]=0
        exhist[-1]=0
        #normalising the histogram to highest point...
        exhist=exhist/(numpy.max(exhist)*1.) 



        width = 0.5
        centers = bBins[:-1]+0.5*(bBins[1:]-bBins[:-1])
        plt.figure()
        #plt.plot(centers, exhist, 'g', linewidth = width)

        from lmfit.models import GaussianModel
        mod = GaussianModel()
        pars = mod.guess(exhist, x=centers)
        out = mod.fit(exhist, pars, x=centers)
        #print (out.fit_report())
        #plt.hist(params, bins = binNo)
        plt.bar(centers,exhist,width=centers[1]-centers[0]) #plot 'real' histogram don't remake it...
        plt.plot(centers,out.best_fit)
        plt.savefig('bootstrap_results/'+param_name+'_boot_hist.pdf')

        #sys.exit(100)
        inty = out.params['amplitude'].value
        av = out.params['center'].value
        sigma = out.params['sigma'].value
        fwhm = out.params['fwhm'].value
        return inty,av,sigma,fwhm


    def plot_results(self,ax, result,BigT, dat, color='none'):
        print((fit_report(result)))
        start = 0.0
        stop = 1.0
        num = len(BigT)
        cm_subsection = numpy.linspace(start,stop,num)
        if color == 'none':
            colors = [cm.jet(x) for x in cm_subsection]
        else:
            colors = []
            for x in range(num):
                colors.append(color)
        fit_vals = self.exchange_func_rest(result.params,  self.BigT_grid, self.q2_grid, BigT)
       
        # for i in range(4):
        for i in range(len(BigT)):
            ax.plot(BigT[i]*(self.gzlvl1*self.delta*self.gyH)**2.0, dat[i,:]/result.params['pre_exp'].value,'x', color = colors[i])
            ax.plot(BigT[i]*(self.gzlvl1*self.delta*self.gyH)**2.0, fit_vals[i,:]/result.params['pre_exp'].value, color = colors[i], label = 'BigT = %f s' %(BigT[i]))
        i = 0.8

        covar_errs = self.getCovar_err(self.param_List, result)

        if(self.mode=='mono'):
            ax.text(0.5, i, r'D$_{out}$ = %1.2e $\pm$ %1.2e cm$^2$ s$^{-1}$ (%s)' %(result.params['D1'].value, covar_errs['D1'], result.params['D1'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'r$_{1,in}$ = %3.2f $\pm$ %3.2f s$^{-1}$ (%s)' %(result.params['r1'].value, covar_errs['r1'], result.params['r1'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03

        elif(self.mode=='bi'):
            #self.param_List = ['D1','D2','pa','r1_1','r1_2','pre_exp']
            ax.text(0.5, i, r'D$_{in}$ = %1.2e $\pm$ %1.2e cm$^2$ s$^{-1}$ (%s)' %(result.params['D2'].value, covar_errs['D2'], result.params['D2'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'D$_{out}$ = %1.2e $\pm$ %1.2e cm$^2$ s$^{-1}$ (%s)' %(result.params['D1'].value, covar_errs['D1'], result.params['D1'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'p$_{out}$ = %0.3e $\pm$ %0.3e (%s)' %(result.params['pa'].value, covar_errs['pa'], result.params['pa'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'r$_{1,in}$ = %3.2f $\pm$ %3.2f s$^{-1}$ (%s)' %(result.params['r1_2'].value, covar_errs['r1_2'], result.params['r1_2'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'r$_{1,out}$ = %3.2f $\pm$ %3.2f s$^{-1}$ (%s)' %(result.params['r1_1'].value, covar_errs['r1_1'], result.params['r1_1'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'I$_0$ = %1.2e $\pm$ %1.2e(%s)' %(result.params['pre_exp'].value, covar_errs['pre_exp'], result.params['pre_exp'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.set_yscale('log')


        else:
            ax.text(0.5, i, r'D$_{in}$ = %1.2e $\pm$ %1.2e cm$^2$ s$^{-1}$ (%s)' %(result.params['D2'].value, covar_errs['D2'], result.params['D2'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'D$_{out}$ = %1.2e $\pm$ %1.2e cm$^2$ s$^{-1}$ (%s)' %(result.params['D1'].value, covar_errs['D1'], result.params['D1'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'p$_{out}$ = %0.3f $\pm$ %0.3f (%s)' %(result.params['pa'].value, covar_errs['pa'], result.params['pa'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'k$_{ex}$ = %4.2f $\pm$ %4.2f s$^{-1}$ (%s)' %(result.params['kex'].value, covar_errs['kex'], result.params['kex'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'r$_{1,in}$ = %3.2f $\pm$ %3.2f s$^{-1}$ (%s)' %(result.params['r1_2'].value, covar_errs['r1_2'], result.params['r1_2'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'r$_{1,out}$ = %3.2f $\pm$ %3.2f s$^{-1}$ (%s)' %(result.params['r1_1'].value, covar_errs['r1_1'], result.params['r1_1'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'I$_0$ = %1.2e $\pm$ %1.2e(%s)' %(result.params['pre_exp'].value, covar_errs['pre_exp'], result.params['pre_exp'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.text(0.5, i, r'radius = %4.2f $\pm$ %4.2f $\mu$m (%s)' %(result.params['rad'].value*1e4, covar_errs['rad']*1e4, result.params['rad'].vary), fontsize = 8,transform=ax.transAxes); i -= 0.03
            ax.set_yscale('log')

        chi = (numpy.sum(self.chiFunc_rest(result.params, self.dat, self.BigT_grid, self.q2_grid, BigT)**2.0))
        ax.text(0.5,i, r'$\chi^2 = %f$' %(chi), fontsize = 8,transform=ax.transAxes)


        ax.set_ylabel('Intensity')
        ax.set_xlabel(r'b = $\Delta(\gamma\delta G)^2$ / cm$^{-2}$ s')
        


    def plot_results_boot(self, boot_av, boot_std ,BigT, dat, plot='none'):
        start = 0.0
        stop = 1.0
        num = len(BigT)
        cm_subsection = numpy.linspace(start,stop,num)
        colors = [cm.jet(x) for x in cm_subsection]

        fit_vals = self.exchange_func_rest(boot_av, self.BigT_grid, self.q2_grid, BigT)
        plt.figure()
        for i in range(len(BigT)):
            plt.plot(BigT[i]*(self.gzlvl1*self.delta*self.gyH)**2.0, dat[i,:],'o', color = colors[i])
            plt.plot(BigT[i]*(self.gzlvl1*self.delta*self.gyH)**2.0, fit_vals[i,:], color = colors[i])
        i = 0.8

        if(self.mode=='mono'):
            plt.figtext(0.5, i, r'D$_{in}$ = %1.2e $\pm$ %1.2e cm$^2$ s$^{-1}$ ' %(boot_av['D1'], boot_std['D1']), fontsize = 8); i -= 0.03
            plt.figtext(0.5, i, r'r$_{1,in}$ = %3.2f $\pm$ %3.2f s$^{-1}$ ' %(boot_av['r1'], boot_std['r1']), fontsize = 8); i -= 0.03
            plt.figtext(0.5, i, r'I$_0$ = %1.2e $\pm$ %1.2e' %(boot_av['pre_exp'], boot_std['pre_exp']), fontsize = 8); i -= 0.03

        else:

            plt.figtext(0.5, i, r'D$_{in}$ = %1.2e $\pm$ %1.2e cm$^2$ s$^{-1}$ ' %(boot_av['D2'], boot_std['D2']), fontsize = 8); i -= 0.03
            plt.figtext(0.5, i, r'D$_{out}$ = %1.2e $\pm$ %1.2e cm$^2$ s$^{-1}$ ' %(boot_av['D1'], boot_std['D1']), fontsize = 8); i -= 0.03
            plt.figtext(0.5, i, r'p$_{out}$ = %0.3f $\pm$ %0.3f ' %(boot_av['pa'], boot_std['pa']), fontsize = 8); i -= 0.03
            plt.figtext(0.5, i, r'k$_{ex}$ = %4.2f $\pm$ %4.2f s$^{-1}$ ' %(boot_av['kex'], boot_std['kex']), fontsize = 8); i -= 0.03
            plt.figtext(0.5, i, r'r$_{1,in}$ = %3.2f $\pm$ %3.2f s$^{-1}$ ' %(boot_av['r1_2'], boot_std['r1_2']), fontsize = 8); i -= 0.03
            plt.figtext(0.5, i, r'r$_{1,out}$ = %3.2f $\pm$ %3.2f s$^{-1}$ ' %(boot_av['r1_1'], boot_std['r1_1']), fontsize = 8); i -= 0.03
            plt.figtext(0.5, i, r'I$_0$ = %1.2e $\pm$ %1.2e' %(boot_av['pre_exp'], boot_std['pre_exp']), fontsize = 8); i -= 0.03
            plt.figtext(0.5, i, r'radius = %4.2f $\pm$ %4.2f $\mu$m' %(boot_av['rad']*1e4, boot_std['rad']*1e4), fontsize = 8); i -= 0.03
            plt.yscale('log')

        boot_chi = (numpy.sum(self.chiFunc_rest(boot_av, self.dat, self.BigT_grid, self.q2_grid, BigT)**2.0))

        plt.figtext(0.5,i, r'$\chi^2 = %f$' %(boot_chi), fontsize = 8)

        

        plt.ylabel('Intensity')
        plt.xlabel(r'b = $\Delta(\gamma\delta G)^2$ / cm$^{-2}$ s')
        if plot != 'none':
            plt.savefig(plot)

    def getCovar_err(self, param_List, result):
        errs = {}
        for item in param_List:
            errs[item] = 0.0
        i = 0
        for item in result.var_names:
            if result.errorbars ==True:
                errs[item] = numpy.sqrt(result.covar[i,i])
            i+=1
        return errs


    def write_results(self, outfile, results):
        stuff = fit_report(results)
        with open(outfile, 'w') as outy:
            outy.write(stuff)

    def write_results_boot(self, outfile, boot_av, boot_std, param_List):
        with open(outfile, 'w') as outy:
            outy.write('bootstrap results\n')
            for item in param_List:
                outy.write('%s average = %s +/- %s\n' %(item, boot_av[item], boot_std[item]))


    def boot_run(self, booty_init,i):
        from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit, fit_report

        print('running boot ', i+1)
        self.fake_var = numpy.zeros_like(self.BigT)
        self.fake_data = numpy.zeros_like(self.dat)

        k = numpy.random.randint(0,len(self.BigT),len(self.BigT)*len(self.gzlvl1))
        l = numpy.random.randint(0, len(self.gzlvl1),len(self.BigT)*len(self.gzlvl1))

        # data with potential replacement for bootstrap
        self.fake_BigT = (self.BigT_grid[k,l]).reshape((self.dat.shape[0],self.dat.shape[1]))
        self.fakeq = (self.q2_grid[k,l]).reshape((self.dat.shape[0],self.dat.shape[1]))
        self.fake_data = (self.dat[k,l]).reshape((self.dat.shape[0],self.dat.shape[1]))


        min_boot = Minimizer(self.chiFunc_rest, booty_init.params, (self.fake_data, self.fake_BigT, self.fakeq, self.BigT), {'boot':True})
        res_boot = min_boot.minimize()
        return res_boot

    def pathExists(self, folder): #create 'out' and 'figs' folders'
        if (os.path.exists(folder)==0):
            os.system('mkdir ' +folder)

    def unpack_boot(self,boot_results):
        if(self.mode=='mono'):

            D1 = []
            r1 = []
            pre_exp = []

            for params in boot_results:
                D1.append(params['D1'].value)
                r1.append(params['r1'].value)
                pre_exp.append(params['pre_exp'].value)


            D1 = numpy.array(D1)
            r1 = numpy.array(r1)
            pre_exp = numpy.array(pre_exp)

            vals_all = {}

            vals_all['D1'] = D1
            vals_all['r1'] = r1
            vals_all['pre_exp'] = pre_exp

            vals_av = {}
            vals_std = {}
            vals_av['D1'] = numpy.average(D1)
            vals_av['r1'] = numpy.average(r1)
            vals_av['pre_exp'] = numpy.average(pre_exp)
            
            vals_std['D1'] = numpy.std(D1)
            vals_std['r1'] = numpy.std(r1)
            vals_std['pre_exp'] = numpy.std(pre_exp)
            return vals_av, vals_std, vals_all


        else:


            D1 = []
            D2 = []
            r1_1 = []
            r1_2 = []
            pa = []
            kex = []
            rad = []
            pre_exp = []

            perm = []
            cell_dens = []

            for params in boot_results:
                D1.append(params['D1'].value)
                D2.append(params['D2'].value)
                r1_1.append(params['r1_1'].value)
                r1_2.append(params['r1_2'].value)
                pa.append(params['pa'].value)
                kex.append(params['kex'].value)
                rad.append(params['rad'].value)
                pre_exp.append(params['pre_exp'].value)

                perm.append(params['kex'].value*params['pa'].value*params['rad'].value/3.0)
                cell_dens.append((1-params['pa'].value)*3/(4*numpy.pi*(params['rad'].value)**3.0))

            D1 = numpy.array(D1)
            D2 = numpy.array(D2)
            r1_1 = numpy.array(r1_1)
            r1_2 = numpy.array(r1_2)
            pa = numpy.array(pa)
            kex = numpy.array(kex)
            rad = numpy.array(rad)
            pre_exp = numpy.array(pre_exp)

            perm = numpy.array(perm)
            cell_dens = numpy.array(cell_dens)

            vals_all = {}

            vals_all['D1'] = D1
            vals_all['D2'] = D2
            vals_all['r1_1'] = r1_1
            vals_all['r1_2'] = r1_2
            vals_all['pa'] = pa
            vals_all['kex'] = kex
            vals_all['rad'] = rad
            vals_all['pre_exp'] = pre_exp


            vals_all['perm'] = perm
            vals_all['cell_dens'] = cell_dens


            vals_av = {}
            vals_std = {}
            vals_av['D1'] = numpy.average(D1)
            vals_av['D2'] = numpy.average(D2)
            vals_av['r1_1'] = numpy.average(r1_1)
            vals_av['r1_2'] = numpy.average(r1_2)
            vals_av['pa'] = numpy.average(pa)
            vals_av['kex'] = numpy.average(kex)
            vals_av['rad'] = numpy.average(rad)
            vals_av['pre_exp'] = numpy.average(pre_exp)
            
            vals_av['perm'] = numpy.average(perm)
            vals_av['cell_dens'] = numpy.average(cell_dens)

            vals_std['D1'] = numpy.std(D1)
            vals_std['D2'] = numpy.std(D2)
            vals_std['r1_1'] = numpy.std(r1_1)
            vals_std['r1_2'] = numpy.std(r1_2)
            vals_std['pa'] = numpy.std(pa)
            vals_std['kex'] = numpy.std(kex)
            vals_std['rad'] = numpy.std(rad)
            vals_std['pre_exp'] = numpy.std(pre_exp)
            
            vals_std['perm'] = numpy.std(perm)
            vals_std['cell_dens'] = numpy.std(cell_dens)
            
            return vals_av, vals_std, vals_all
