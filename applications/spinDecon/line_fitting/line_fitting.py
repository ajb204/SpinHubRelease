#!/usr/bin/python
import numpy
import matplotlib            #import matplotlib
matplotlib.use('WXAgg')      #switch on the wxPython mode
from matplotlib.widgets import Cursor, RectangleSelector
from spinDecon.gui.dialogs.processing.process import path_escape
from spinDecon.processing.nmrpipe_scripts import MakeProj4D, MakeProj3D
import scipy.optimize as opt
import hashlib, pickle
import os, re


def Gauss(x,x0,sigma):
   
    return (1/(sigma*numpy.sqrt(2*numpy.pi)))*numpy.exp(-(x-x0)**2./(2*sigma**2.)) #*1./(numpy.sqrt(2*numpy.pi)*sigma)

def Lorentz(x,x0,Gamma):
    nom = Gamma/2.
    denom = (x-x0)**2.+(Gamma/2.)**2.
    norm = 1/(numpy.pi)
   
    return norm*nom/denom

def PV(amp,xx,x0,Gamma,sigma,nu):
    yvals=(1-nu)*Gauss(xx,x0,sigma)+(nu)*(Lorentz(xx,x0,Gamma))
    return amp*(yvals) #/numpy.max(yvals))


def Overlap_PseudoVoigt(*params):
    length = len(params)
    output = []
    if (length-8) % 3 == 0:
        overlaps = int((length-7)/3.)
        coords=params[0]
        noiseval = params[1]
        nu2 = params[-1]
        nu1 = params[-2]
        sigma_y = params[-3]
        sigma_x= params[-4]
        Gamma_y = params[-5]
        Gamma_x = params[-6]
        for x in range(overlaps):
            amp = params[1+(3*x)]
            x0 = params[2+(3*x)]
            y0 = params[3+(3*x)]

            if len(output) == 0:
                output = PseudoVoigt(coords, 0, amp, x0,y0,Gamma_x, Gamma_y, sigma_x, sigma_y, nu1, nu2)

            else:
                output+=PseudoVoigt(coords, 0, amp, x0,y0,Gamma_x, Gamma_y, sigma_x, sigma_y, nu1, nu2)
        # exit()
    else:
        print(length)
        exit()

    return output*(output>noiseval)

def Overlap_PseudoVoigt_individual(*params):
    length = len(params)
    output = []
    # print(params)
    if (length-2) % 9 == 0:
        overlaps = int((length-1)/9.)
        coords=params[0]
        noiseval = params[1]
        for x in range(overlaps):
            amp = params[2+(9*x)]
            x0 = params[3+(9*x)]
            y0 = params[4+(9*x)]
            Gamma_x = params[5+(9*x)]
            Gamma_y = params[6+(9*x)]
            sigma_x = params[7+(9*x)]
            sigma_y = params[8+(9*x)]
            nu1 = params[9+(9*x)]
            nu2 = params[10+(9*x)]

            if len(output) == 0:
                output = PseudoVoigt(coords, 0, amp, x0,y0,Gamma_x, Gamma_y, sigma_x, sigma_y, nu1, nu2)
            else:
                output+=PseudoVoigt(coords, 0, amp, x0,y0,Gamma_x, Gamma_y, sigma_x, sigma_y, nu1, nu2)
        # exit()
    else:
        print(length)
        exit()

    return output*(output>noiseval)

def PseudoVoigt(coords, noiseval, amp, x0,y0,Gamma_x, Gamma_y, sigma_x, sigma_y, nu1, nu2):
    x,y=coords
    
    yvals=((1-nu1)*Gauss(x,x0,sigma_x)+(nu1)*(Lorentz(x,x0,Gamma_x)))*((1-nu2)*Gauss(y,y0,sigma_y)+(nu2)*(Lorentz(y,y0,Gamma_y)))
    ans = amp*(yvals) #/numpy.max(yvals))
    # print((ans[ans>noiseval]))
    ans = ans*(numpy.fabs(ans)>noiseval)
    
    return ans.ravel()



class Unidec_line_fitting():

    def __init__(self, data, peak_list, peak_list_names, gamma_x, gamma_y, sigma_x, sigma_y, nu1, nu2, frac_thresh, uc0, uc1):
        self.Gamma_x = float(gamma_x)/(2.*float(uc0.ppm_scale()[0]-uc0.ppm_scale()[1]))
        self.Gamma_y = float(gamma_y)/(2.*float(uc1.ppm_scale()[0]-uc1.ppm_scale()[1]))
        self.sigma_x = float(sigma_x)/(2.*float(uc0.ppm_scale()[0]-uc0.ppm_scale()[1]))
        self.sigma_y = float(sigma_y)/(2.*float(uc1.ppm_scale()[0]-uc1.ppm_scale()[1]))

        print(self.Gamma_x, self.Gamma_y, self.sigma_x, self.sigma_y)
        # exit()
        self.nu1 = float(nu1)
        self.nu2 = float(nu2)
        self.data = data
        if len(self.data.shape)==2:
            self.data = numpy.expand_dims(self.data, axis=0)
        self.peak_list = peak_list
        self.peak_list_names = peak_list_names
        self.plotting_fuda_data = {}
        self.plotting_resim_data = {}
        self.overlap_resim_data = {}
        self.fitted_peaks = []
        self.unoverlapped = []
        self.unoverlapped_names = []
        self.fitted_names = []
        self.nu2_array = []
        self.nu1_array = []
        self.sigma_y_array = []
        self.sigma_x_array = []
        self.Gamma_y_array = []
        self.Gamma_x_array = []
        self.overlap_distance = 3.0 ## Define the overlap distance as a multiple
        self.frac_thresh = frac_thresh
        self.thresh = 0. #self.frac_thresh*numpy.max(numpy.abs(self.data))
        self.intensities = {}
        self.uc0 = uc0
        self.uc1 = uc1
        # self.fwhm_x_average = 0




    def finding_overlaps(self):
        
        i = 0

        # self.overlapped_dic = {}
        overlap_array=[]
        for number, peak in enumerate(self.peak_list):
            x_int = int(peak[0])
            y_int = int(peak[1])
            
            # label = peak[4]

            overlapped=False ## We are checking for overlap
            peak_list2 = numpy.delete(self.peak_list, i, axis=0) ## remove current peak from peak list?
            peak_list2_names = numpy.delete(self.peak_list_names, i, axis=0) ## remove current peak from peak list?

            # j=0
            args = numpy.argwhere(numpy.abs(x_int-peak_list2[:,0])< self.fwhm_x*self.overlap_distance) ## Returning the other peaks within the FWHM of the x dimension
            args2 = numpy.squeeze(numpy.argwhere(numpy.abs(y_int-peak_list2[numpy.squeeze(args),1])< self.fwhm_y*self.overlap_distance)) ## Returning the other peaks within the FWHM of both the x and y dimensions (i.e. coords of the closest peaks)

            x_yes = peak_list2[args[args2]] ## returning the actual peaks that are overlapped
            x_yes_names = peak_list2_names[args[args2]]
            if len(x_yes.shape) == 2: ## making the array dimensions consistent
                x_yes = [x_yes]



            for x_y_number, x_y in enumerate(x_yes): ## Iterate through overlapped peaks


                # print(int(args[args2][x_y_number]))
                x_y= x_y[0] ## Remove the first indent of the array
                overlap_array.append([peak, x_y, self.peak_list_names[number], numpy.squeeze(x_yes_names[x_y_number])]) ## Append 'the overlap' to the overlap array
                overlapped = True ## The current peak is overlapped!
            # exit()
            ## If the peak is not overlapped, append it to the unoverlapped array.
            if overlapped == False:
                y_1 = max(0,int(x_int-self.fwhm_x*self.overlap_distance))
                y_2 = min(self.data.shape[1],int(x_int+self.fwhm_x*self.overlap_distance))
                x_1 = max(0,int(y_int-self.fwhm_y*self.overlap_distance))
                x_2 = min(self.data.shape[2],int(y_int+self.fwhm_y*self.overlap_distance))
                ppm_x = self.uc1.ppm(y_int)
                ppm_y = self.uc0.ppm(x_int)
                self.unoverlapped.append([y_1,y_2, x_1,x_2, ppm_x, ppm_y, y_int, x_int, self.peak_list_names[number]])
                self.unoverlapped_names.append(str(self.peak_list_names[number]))
            i+=1

        self.final_overlaps=[[]]
        self.final_overlap_names=[[]]
        x=0
        # print(overlap_array)
        # exit()

        for x in overlap_array: ## Iterating through array of overlaps
            pk1 = x[0] ## First peak in the overlap
            pk2 = x[1] ## Second peak in the overlap
            name1 = str(x[2])
            name2 = str(x[3])
            overlapping = False ## Default position
            for ol in range(len(self.final_overlap_names)): ## Iterate through the islands to see if either peak is there
                if name1 in self.final_overlap_names[ol]: ## Check if first peak is in the array
                    if name2 not in self.final_overlap_names[ol]: ## If the second peak is not in the array, append it
                        self.final_overlaps[ol].append(pk2)
                        self.final_overlap_names[ol].append(name2)
                        overlapping=True
                    else: ## If second peak is already there, we don't want to be in the default position.
                        overlapping=True

                elif name2 in self.final_overlap_names[ol]: ## Check if second peak is in the array but first isnt: if so append to island!
                    self.final_overlaps[ol].append(pk1)
                    self.final_overlap_names[ol].append(name1)
                    overlapping=True

            if overlapping==False: ## If we are still in the default position, begin the island
                self.final_overlaps.append([pk1, pk2])
                self.final_overlap_names.append([name1, name2])

        del(self.final_overlap_names[0])
        del(self.final_overlaps[0])

        self.total_unoverlapped=len(self.unoverlapped)


    def plot_fuda_fit(self, peak, number, ax, canvas):
        ax.cla()

        popt = self.plotting_resim_data[peak][number]
        popt[0] = 0



        x = numpy.arange(len(self.plotting_fuda_data[peak][number][0,:]))
        y = numpy.arange(len(self.plotting_fuda_data[peak][number][:,0]))
        x,y = numpy.meshgrid(x,y)
        print(len(popt))
        if peak in self.overlap_resim_data.keys():
            resim = numpy.reshape(Overlap_PseudoVoigt_individual((x,y),*popt), self.plotting_fuda_data[peak][number][:,:].shape)
            popt2 = self.overlap_resim_data[peak][number]
            popt2[0] = 0

            resim2 = numpy.reshape(Overlap_PseudoVoigt_individual((x,y),*popt2), self.plotting_fuda_data[peak][number][:,:].shape)
            ax.plot_wireframe(x,y,resim2, color='g', zorder=500, alpha=0.5)
        else:
            # print(popt)
            print(popt[1], popt[4], popt[5], popt[8], popt[9])
            resim = numpy.reshape(PseudoVoigt((x,y),*popt), self.plotting_fuda_data[peak][number][:,:].shape)


        ax.plot_wireframe(x,y, self.plotting_fuda_data[peak][number][:,:])
        ax.plot_wireframe(x,y,resim, color='r', zorder=1000)
        ax.plot_wireframe(x,y,numpy.ones_like(x)*self.thresh, color='k', zorder=1000, alpha=0.25)
        # ax.plot_wireframe(x,y,self.plotting_fuda_data[peak][number][:,:]*(self.plotting_fuda_data[peak][number][:,:]>self.thresh), color='g', zorder=1000, alpha=0.5)
        # ax.plot_wireframe(x,y,numpy.ones_like(x)*numpy.max(numpy.abs(self.plotting_fuda_data[peak][number][:,:]))*self.frac_thresh, color='k', zorder=1000, alpha=0.5)
        parameter_string = 'FWHM x = %.2f ppm (av = %.2f)\nFWHM y = %.2f ppm (av = %.2f)' % (self.intensities[peak][1][number], self.fwhm_x_average,self.intensities[peak][2][number], self.fwhm_y_average)
        ax.text2D(0.05, 0.95,parameter_string, transform=ax.transAxes)
        # self.bar_chart(bar_ax)
        canvas.draw()

    def fit_unoverlapped_peak(self, name, ftol=1e-12):
        for x, peak_name in enumerate(self.unoverlapped_names):
            if name == peak_name:
                peak = self.unoverlapped[x]
                # if peak[0] > 0 and peak[2]>0:
                fuda_data = self.data[:,peak[0]:peak[1], peak[2]:peak[3]]
                intensity, fwhm_x, fwhm_y, amp = self.fuda(fuda_data, (numpy.round((peak[1]-peak[0])/2.),numpy.round((peak[3]-peak[2])/2.)), name = str(peak[8]), ftol=ftol)
                self.fitted_peaks.append([peak[5], peak[4]])
                self.fitted_names.append([peak[8]])
                self.intensities[peak[8]] = [intensity, fwhm_x, fwhm_y, amp]
                self.calculate_averages()
                self.print_intensities()
                self.save_results()
                return

    

    def is_peak_overlapped(self, name):
        if name in self.peak_list_names:
            if name in self.unoverlapped_names:
                return False
            else:
                for x, island in enumerate(self.final_overlap_names):
                    if name in island:
                        return x
        else:
            print('Error: Peak not found in main peak list')
            return None

    def fit_overlapped_peaks(self, x, ftol=1e-10):
        peaks, data = self.determine_overlap_area(self.final_overlaps[x])
        print(self.final_overlap_names[x], peaks)
        integ, fwhm_x, fwhm_y, amp = self.overlap_fuda(data, peaks, self.final_overlap_names[x], ftol=ftol)
        for inte in range(integ.shape[1]):
            self.intensities[self.final_overlap_names[x][inte]] = [integ[:,inte],fwhm_x[:,inte],fwhm_y[:,inte], amp[:,inte]]
        self.fitted_names.append(self.final_overlap_names[x])
        self.fitted_peaks.append([self.final_overlaps[x][0][2], self.final_overlaps[x][0][3]])
        self.calculate_averages()

        self.print_intensities()
        self.save_results


    def calculate_averages(self):
        self.fwhm_x_average = 0
        self.fwhm_y_average = 0
        count = 0
        print('calculating')
        print(self.intensities.keys())
        for key in self.intensities.keys():
            print(self.intensities[key])
            self.fwhm_x_average += numpy.mean(self.intensities[key][1])
            self.fwhm_y_average += numpy.mean(self.intensities[key][2])
            count += 1
        self.fwhm_x_average = float(self.fwhm_x_average/float(count))
        self.fwhm_y_average = float(self.fwhm_y_average/float(count))
        

    def prelim_fuda_thread(self, data_coords, name):

        fuda_data = self.data[:,data_coords[0]:data_coords[1], data_coords[2]:data_coords[3]]

        intensity, fwhm_x, fwhm_y, amp = self.fuda(fuda_data, numpy.unravel_index(fuda_data[0,:,:].argmax(), fuda_data[0,:,:].shape), name)
        # self.scatter_data.append(intensity)


        # self.plot_scatters(self.axes_scatter)
        av_vals = numpy.array(self.resim_data[-1])
        resim_x = PV(1, numpy.linspace(0,self.fuda_data[0].shape[0], 200), self.fuda_data[0].shape[0]/2., av_vals[4], av_vals[6], av_vals[-2])[:101] # Calculate the fwhm of the peak and then only take the first half (we have two minima)
        resim_y = PV(1, numpy.linspace(0,self.fuda_data[0].shape[1], 200), self.fuda_data[0].shape[1]/2., av_vals[5], av_vals[7], av_vals[-1])[:101] # Calculate the fwhm of the peak and then only take the first half (we have two minima)
        frac_hwhm_x = 100-numpy.argmin(numpy.fabs(resim_x-(0.5*numpy.max(numpy.fabs(resim_x))))) # how far away are we from the centre of the PV function in simulated PV units?
        frac_hwhm_y = 100-numpy.argmin(numpy.fabs(resim_y-(0.5*numpy.max(numpy.fabs(resim_y))))) # how far away are we from the centre of the PV function in simulated PV units?
        hwhm_x=numpy.abs(fuda_data[0].shape[0]*frac_hwhm_x/200.) # convert between simulated PV units and actual data units
        hwhm_y=numpy.abs(fuda_data[0].shape[1]*frac_hwhm_y/200.) # convert between simulated PV units and actual data units
        # self.fwhm_x = hwhm_y*2
        # self.fwhm_y = hwhm_x*2
        
        self.fwhm_x = 32.19
        self.fwhm_y = 17.76

        # self.fwhm_x = fwhm_x
        # self.fwhm_y = fwhm_y
        # print(self.fwhm_x, self.fwhm_y)
        # exit()
        self.intensities[name] = [intensity, fwhm_x, fwhm_y, amp]
        self.calculate_averages()
        return intensity

    def determine_overlap_area(self, peaks):
        lowest_y = 1e10
        highest_y = -1
        lowest_x = 1e10
        highest_x = -1
        mask = numpy.ones_like(self.data)
        peaks2=[]
        # names = []
        for numbers, peak in enumerate(peaks):
            # names.append(peak[4])
            y = peak[1]
            x = peak[0]
            x_1 = int(x-self.fwhm_x*self.overlap_distance*1./2.)
            x_2 = int(x+self.fwhm_x*self.overlap_distance*1./2.)
            y_1 = int(y-self.fwhm_y*self.overlap_distance*1./2.)
            y_2 = int(y+self.fwhm_y*self.overlap_distance*1./2.)
            if y_1 < lowest_y:
                lowest_y = y_1
            if x_1 < lowest_x:
                lowest_x = x_1
            if y_2 > highest_y:
                highest_y = y_2
            if x_2 > highest_x:
                highest_x = x_2
            mask[:,x_1:x_2, y_1:y_2] = 0

        for peak in peaks:
            y = peak[1]
            x = peak[0]
            peaks2.append([y-lowest_y, x-lowest_x])

        masked = numpy.ma.masked_array(self.data, mask=mask, fill_value=0.0)[:,lowest_x:highest_x, lowest_y:highest_y]
        return peaks2, masked.filled()

    def print_intensities(self):
        filename = 'out/intensity_fitting.out'
        outy = open(filename, 'w')

        for key in self.intensities.keys():
            outy.write("Peak Name: "+str(key)+"\n")

            for x in self.intensities[key][0]:
                outy.write("%f\n" % x)
            outy.write("\n")
        outy.close()

    # def overlap_fuda_fit(self):

    def save_overlap_values(self, popt):
        if len(self.nu2_array)>0:
            if numpy.fabs((popt[-1]-self.nu2)/self.nu2) < 3.0:
                self.nu2_array.append(popt[-1])
                print("saving nu2:", popt[-1])

            if numpy.fabs((popt[-2]-self.nu1)/self.nu1) < 3.0:
                self.nu1_array.append(popt[-2])
                print("saving nu1:", popt[-2])

            if numpy.fabs((popt[-3]-self.sigma_y)/self.sigma_y) < 3.0:
                self.sigma_y_array.append(popt[-3])
                print("saving sigma_y:", popt[-3])

            if numpy.fabs((popt[-4]-self.sigma_x)/self.sigma_x) < 3.0:
                self.sigma_x_array.append(popt[-4])
                print("saving sigma_x:", popt[-4])

            if numpy.fabs((popt[-5]-self.Gamma_y)/self.Gamma_y) < 3.0:
                self.Gamma_y_array.append(popt[-5])
                print("saving gamma_y:", popt[-5])

            if numpy.fabs((popt[-6]-self.Gamma_x)/self.Gamma_x) < 3.0:
                self.Gamma_x_array.append(popt[-6])
                print("saving gamma_x:", popt[-6])

        else:
            self.nu2_array.append(popt[-1])
            self.nu1_array.append(popt[-2])
            self.sigma_y_array.append(popt[-3])
            self.sigma_x_array.append(popt[-4])
            self.Gamma_y_array.append(popt[-5])
            self.Gamma_x_array.append(popt[-6])


        self.nu2 = numpy.average(self.nu2_array)
        self.nu1 = numpy.average(self.nu1_array)
        self.sigma_y = numpy.average(self.sigma_y_array)
        self.sigma_x = numpy.average(self.sigma_x_array)
        self.Gamma_y = numpy.average(self.Gamma_y_array)
        self.Gamma_x = numpy.average(self.Gamma_x_array)

    def save_overlap_values_individual(self, params):
        length = len(params)
        if (length-1) % 9 == 0:
            overlaps = int((length-1)/9.)
            
            for x in range(overlaps):
                Gamma_x = params[4+(3*x)]
                Gamma_y = params[5+(3*x)]
                sigma_x = params[6+(3*x)]
                sigma_y = params[7+(3*x)]
                nu1 = params[8+(3*x)]
                nu2 = params[9+(3*x)]
            

            if len(self.nu2_array)>0:
                if numpy.fabs((nu2-self.nu2)/self.nu2) < 3.0:
                    self.nu2_array.append(nu2)
                    print("saving nu2:", nu2)
                    
                if numpy.fabs((nu1-self.nu1)/self.nu1) < 3.0:
                    self.nu1_array.append(nu1)
                    print("saving nu1:", nu1)

                if numpy.fabs((sigma_y-self.sigma_y)/self.sigma_y) < 3.0:
                    self.sigma_y_array.append(sigma_y)
                    print("saving sigma_y:", sigma_y)

                if numpy.fabs((sigma_x-self.sigma_x)/self.sigma_x) < 3.0:
                    self.sigma_x_array.append(sigma_x)
                    print("saving sigma_x:", sigma_x)

                if numpy.fabs((Gamma_y-self.Gamma_y)/self.Gamma_y) < 3.0:
                    self.Gamma_y_array.append(Gamma_y)
                    print("saving gamma_y:", Gamma_y)

                if numpy.fabs((Gamma_x-self.Gamma_x)/self.Gamma_x) < 3.0:
                    self.Gamma_x_array.append(Gamma_x)
                    print("saving gamma_x:", Gamma_x)

            else:
                self.nu2_array.append(nu2)
                self.nu1_array.append(nu1)
                self.sigma_y_array.append(sigma_y)
                self.sigma_x_array.append(sigma_x)
                self.Gamma_y_array.append(Gamma_y)
                self.Gamma_x_array.append(Gamma_x)


        self.nu2 = numpy.average(self.nu2_array)
        self.nu1 = numpy.average(self.nu1_array)
        self.sigma_y = numpy.average(self.sigma_y_array)
        self.sigma_x = numpy.average(self.sigma_x_array)
        self.Gamma_y = numpy.average(self.Gamma_y_array)
        self.Gamma_x = numpy.average(self.Gamma_x_array)
        

    def save_fuda_values(self, popt):
        if len(self.nu2_array)>0:
            if numpy.fabs((popt[-1]-self.nu2)/self.nu2) < 3.0:
                self.nu2_array.append(popt[-1])
                print("saving nu2:", popt[-1])
                
            if numpy.fabs((popt[-2]-self.nu1)/self.nu1) < 3.0:
                self.nu1_array.append(popt[-2])
                print("saving nu1:", popt[-2])

            if numpy.fabs((popt[-3]-self.sigma_y)/self.sigma_y) < 3.0:
                self.sigma_y_array.append(popt[-3])
                print("saving sigma_y:", popt[-3])

            if numpy.fabs((popt[-4]-self.sigma_x)/self.sigma_x) < 3.0:
                self.sigma_x_array.append(popt[-4])
                print("saving sigma_x:", popt[-4])

            if numpy.fabs((popt[-5]-self.Gamma_y)/self.Gamma_y) < 3.0:
                self.Gamma_y_array.append(popt[-5])
                print("saving gamma_y:", popt[-5])

            if numpy.fabs((popt[-6]-self.Gamma_x)/self.Gamma_x) < 3.0:
                self.Gamma_x_array.append(popt[-6])
                print("saving gamma_x:", popt[-6])

        else:
            self.nu2_array.append(popt[-1])
            print("saving nu2:", popt[-1])
            self.nu1_array.append(popt[-2])
            print("saving nu1:", popt[-2])
            self.sigma_y_array.append(popt[-3])
            print("saving sigma_y:", popt[-3])
            self.sigma_x_array.append(popt[-4])
            print("saving sigma_x:", popt[-4])
            self.Gamma_y_array.append(popt[-5])
            print("saving gamma_y:", popt[-5])
            self.Gamma_x_array.append(popt[-6])
            print("saving gamma_x:", popt[-6])


        self.nu2 = numpy.average(self.nu2_array)
        self.nu1 = numpy.average(self.nu1_array)
        self.sigma_y = numpy.average(self.sigma_y_array)
        self.sigma_x = numpy.average(self.sigma_x_array)
        self.Gamma_y = numpy.average(self.Gamma_y_array)
        self.Gamma_x = numpy.average(self.Gamma_x_array)

    def overlap_fuda(self, full_data, peaks,names, ftol):
        if len(full_data.shape)==2:
            full_data = numpy.expand_dims(full_data, axis=0)
        x = numpy.arange(len(full_data[0,0,:]))
        y = numpy.arange(len(full_data[0,:,0]))
        x,y = numpy.meshgrid(x,y)

        self.fuda_data = full_data
        self.resim_data = []

        intensities = []
        fwhm_x_finals = []
        fwhm_y_finals = []
        amp_finals=[]
        plotting_resim_data_row = []
        plotting_fuda_data_row = []

        fudas = {}
        resims = {}

        print(names)
        # exit()

        for x2 in range(len(peaks)):
            name = names[x2]
            resims[name] = []
            fudas[name] = []

        if len(self.fuda_data.shape) == 2:
            self.fuda_data = [self.fuda_data]

        for number,dat in enumerate(self.fuda_data):
            initial_guess = []
            bounds_neg = []
            bounds_pos = []
            integrals = []
            fwhm_x_array=[]
            fwhm_y_array=[]
            amps = []

            self.thresh = float(numpy.max(numpy.fabs(dat))*self.frac_thresh)

            if number == 0:
                initial_guess.append(self.thresh)
                bounds_neg.append(self.thresh-1e-10)
                bounds_pos.append(self.thresh+1e-10)
                for x2 in range(len(peaks)):
                    amp_factor  = (numpy.reshape(PseudoVoigt((x,y), self.thresh, 1., int(peaks[x2][0]), int(peaks[x2][1]), self.Gamma_x, self.Gamma_y, self.sigma_x, self.sigma_y, self.nu1, self.nu2), dat.shape)[int(peaks[x2][1]), int(peaks[x2][0])])
                    initial_guess.append(dat[int(peaks[x2][1]), int(peaks[x2][0])]/amp_factor)
                    # initial_guess.append(dat[int(peaks[x2][1]), int(peaks[x2][0])])
                    bounds_neg.append(-numpy.inf)
                    bounds_pos.append(numpy.inf)
                    initial_guess.append(peaks[x2][0])
                    bounds_neg.append(-numpy.inf)
                    bounds_pos.append(numpy.inf)
                    initial_guess.append(peaks[x2][1])
                    bounds_neg.append(-numpy.inf)
                    bounds_pos.append(numpy.inf)
                self.default_values_gauss(initial_guess, bounds_neg, bounds_pos)

                initial_guess = list(initial_guess)
                popt, pcov = opt.curve_fit(Overlap_PseudoVoigt, (x,y), (dat*(numpy.abs(dat)>self.thresh)).ravel(), p0 = initial_guess, maxfev=100000,ftol=ftol,gtol=1e-10,xtol=1e-12, verbose = 1, bounds = (bounds_neg, bounds_pos))

                self.save_overlap_values(popt)
                initial_guess = []
                bounds_neg = []
                bounds_pos = []
                initial_guess.append(self.thresh)
                bounds_neg.append(self.thresh-1e-10)
                bounds_pos.append(self.thresh+1e-10)
                for x2 in range(len(peaks)):
                    amp_factor  = (numpy.reshape(PseudoVoigt((x,y), 0., 1., int(peaks[x2][0]), int(peaks[x2][1]), self.Gamma_x, self.Gamma_y, self.sigma_x, self.sigma_y, self.nu1, self.nu2), dat.shape)[int(peaks[x2][1]), int(peaks[x2][0])])
                    initial_guess.append(dat[int(peaks[x2][1]), int(peaks[x2][0])]/amp_factor)
                    bounds_neg.append(-numpy.inf)
                    bounds_pos.append(numpy.inf)
                    initial_guess.append(peaks[x2][0])
                    bounds_neg.append(-numpy.inf)
                    bounds_pos.append(numpy.inf)
                    initial_guess.append(peaks[x2][1])
                    bounds_neg.append(-numpy.inf)
                    bounds_pos.append(numpy.inf)
                    self.default_values_gauss(initial_guess, bounds_neg, bounds_pos)

                print('Relaxing peakshape for each peak')
                initial_guess = list(initial_guess)
                popt, pcov = opt.curve_fit(Overlap_PseudoVoigt_individual, (x,y), (dat*(numpy.abs(dat)>self.thresh)).ravel(), p0 = initial_guess, maxfev=1000,ftol=ftol,gtol=1e-10,xtol=1e-10, verbose = 1, bounds = (bounds_neg, bounds_pos))

                self.save_overlap_values_individual(popt)

            else:

                for x2 in range(len(peaks)):
                    initial_guess.append(self.thresh)
                    bounds_neg.append(self.thresh-1e-10)
                    bounds_pos.append(self.thresh+1e-10)
                    amp_factor  = (numpy.reshape(PseudoVoigt((x,y), 0., 1., int(peaks[x2][0]), int(peaks[x2][1]), self.Gamma_x, self.Gamma_y, self.sigma_x, self.sigma_y, self.nu1, self.nu2), dat.shape)[int(peaks[x2][1]), int(peaks[x2][0])])
                    initial_guess.append(dat[int(peaks[x2][1]), int(peaks[x2][0])]/amp_factor)
                    bounds_neg.append(-numpy.inf)
                    bounds_pos.append(numpy.inf)
                    initial_guess.append(peaks[x2][0])
                    bounds_neg.append(peaks[x2][0]-2.)
                    bounds_pos.append(peaks[x2][0]+2.)
                    initial_guess.append(peaks[x2][1])
                    bounds_neg.append(peaks[x2][1]-2.)
                    bounds_pos.append(peaks[x2][1]+2.)
                    initial_guess.append(self.Gamma_x)
                    bounds_neg.append(self.Gamma_x-0.25*numpy.abs(self.Gamma_x))
                    bounds_pos.append(self.Gamma_x+0.25*numpy.abs(self.Gamma_x))
                    initial_guess.append(self.Gamma_y)
                    bounds_neg.append(self.Gamma_y-0.25*numpy.abs(self.Gamma_y))
                    bounds_pos.append(self.Gamma_y+0.25*numpy.abs(self.Gamma_y))
                    initial_guess.append(self.sigma_x)
                    bounds_neg.append(self.sigma_x-0.25*numpy.abs(self.sigma_x))
                    bounds_pos.append(self.sigma_x+0.25*numpy.abs(self.sigma_x))
                    initial_guess.append(self.sigma_y)
                    bounds_neg.append(self.sigma_y-0.25*numpy.abs(self.sigma_y))
                    bounds_pos.append(self.sigma_y+0.25*numpy.abs(self.sigma_y))
                    initial_guess.append(self.nu1)
                    bounds_neg.append(0)
                    bounds_pos.append(1)
                    initial_guess.append(self.nu2)
                    bounds_neg.append(0)
                    bounds_pos.append(1)

                initial_guess = list(initial_guess)
                popt, pcov = opt.curve_fit(Overlap_PseudoVoigt_individual, (x,y), (dat*(numpy.abs(dat)>self.thresh)).ravel(), p0 = initial_guess, maxfev=100000,ftol=ftol,gtol=1e-10,xtol=1e-12, bounds = (bounds_neg, bounds_pos))

            self.resim_data.append(popt)

            plotting_resim_data_row.append(popt)
            plotting_fuda_data_row.append(dat)

            initial_guess2 = []
            bounds_neg2 = []
            bounds_pos2 = []

            for x2 in range(len(peaks)):
                name = names[x2]
                amp = popt[1+(9*x2)]
                x0 = popt[2+(9*x2)]
                y0 = popt[3+(9*x2)]
                gamma_x = popt[4+(9*x2)]
                gamma_y = popt[5+(9*x2)]
                sigma_x = popt[6+(9*x2)]
                sigma_y = popt[7+(9*x2)]
                nu1 = popt[8+(9*x2)]
                nu2 = popt[9+(9*x2)]

                resim_each = (0, amp, x0, y0, gamma_x, gamma_y, sigma_x, sigma_y, nu1, nu2)
                integrals.append(numpy.sum(PseudoVoigt((x,y),*resim_each)))
                # integrals.append(1)
                resims[name].append(resim_each)
                fudas[name].append(dat)
                amps.append(amp)


                initial_guess2.append(amp) #intensity
                bounds_neg2.append(-numpy.inf)
                bounds_pos2.append(numpy.inf)
                initial_guess2.append(x0) # x0
                bounds_neg2.append(-numpy.inf)
                bounds_pos2.append(numpy.inf)
                initial_guess2.append(y0) # y0
                bounds_neg2.append(-numpy.inf)
                bounds_pos2.append(numpy.inf)
                initial_guess2.append(gamma_x) #gamma x
                bounds_neg2.append(-numpy.inf)
                bounds_pos2.append(numpy.inf)
                initial_guess2.append(gamma_y) #gamma y
                bounds_neg2.append(-numpy.inf)
                bounds_pos2.append(numpy.inf)
                initial_guess2.append(sigma_x) #sigma x
                bounds_neg2.append(-numpy.inf)
                bounds_pos2.append(numpy.inf)
                initial_guess2.append(sigma_y) #sigma y
                bounds_neg2.append(-numpy.inf)
                bounds_pos2.append(numpy.inf)
                initial_guess2.append(nu1) #nu1
                bounds_neg2.append(0)
                bounds_pos2.append(1)
                initial_guess2.append(nu2) #nu2
                bounds_neg2.append(0)
                bounds_pos2.append(1)

                
                

                
            fwhm_x_finals.append(fwhm_x_array)
            fwhm_y_finals.append(fwhm_y_array)
            amp_finals.append(amps)
            intensities.append(integrals)

        for key in resims.keys():
            self.plotting_resim_data[key] = numpy.array(resims[key])
            self.plotting_fuda_data[key] = numpy.array(fudas[key])
            self.overlap_resim_data[key] = numpy.array(self.resim_data)


            resim_x = PV(1, numpy.linspace(0,self.fuda_data.shape[0], 200), self.fuda_data.shape[0]/2., self.plotting_resim_data[key][0][4], self.plotting_resim_data[key][0][6], self.plotting_resim_data[key][0][-2])[:101]
            resim_y = PV(1, numpy.linspace(0,self.fuda_data.shape[1], 200), self.fuda_data.shape[1]/2., self.plotting_resim_data[key][0][5], self.plotting_resim_data[key][0][7], self.plotting_resim_data[key][0][-1])[:101]
            frac_hwhm_x = 100-numpy.argmin(numpy.fabs(resim_x-(0.5*numpy.max(numpy.fabs(resim_x))))) # how far away are we from the centre of the PV function in simulated PV units?
            frac_hwhm_y = 100-numpy.argmin(numpy.fabs(resim_y-(0.5*numpy.max(numpy.fabs(resim_y))))) # how far away are we from the centre of the PV function in simulated PV units?
            hwhm_x=numpy.abs(self.fuda_data[0].shape[0]*frac_hwhm_x/200.) # convert between simulated PV units and actual data units
            hwhm_y=numpy.abs(self.fuda_data[0].shape[1]*frac_hwhm_y/200.) # convert between simulated PV units and actual data units
            fwhm_x = hwhm_y*2
            fwhm_y = hwhm_x*2

            # self.fwhm_x = fwhm_x
            # self.fwhm_y = fwhm_y

            fwhm_x_array.append(-self.uc1.ppm(fwhm_x)+ self.uc1.ppm(0))
            fwhm_y_array.append(-self.uc0.ppm(fwhm_y)+ self.uc0.ppm(0))

        return numpy.array(intensities), numpy.array(fwhm_x_finals), numpy.array(fwhm_y_finals), numpy.array(amp_finals)


    def save_results(self):
        
        current_correlate_hash = self.hash_correlate()
        
        save_array = [current_correlate_hash, self.plotting_resim_data, self.plotting_fuda_data, self.intensities, self.final_overlaps, self.final_overlap_names, self.fwhm_x, self.fwhm_y]
        pickle.dump(save_array, open('out/fuda.save', 'wb'))

        

    def hash_correlate(self):
        BUF_SIZE = 65536  # lets read stuff in 64kb chunks!

        md5 = hashlib.md5()
        # sha1 = hashlib.sha1()
        try:
            with open('out/correlate.3', 'rb') as f:
                while True:
                    data = f.read(BUF_SIZE)
                    if not data:
                        break
                    md5.update(data)
                    # sha1.update(data)
            
        except:
            return 0
        return str(md5.hexdigest())



    def load_results(self):
        if os.path.exists('out/fuda.save'):
            save_array = pickle.load(open('out/fuda.save', 'rb'))
            saved_hash = save_array[0]
            if saved_hash != self.hash_correlate():
                print('We have a save file but the correlate.3 file does not match up: redo!')
                return False
                
            else:
                try:
                    self.plotting_resim_data = save_array[1]
                    self.plotting_fuda_data = save_array[2]
                    self.intensities = save_array[3]
                    self.final_overlaps = save_array[4]
                    self.final_overlap_names = save_array[5]
                    self.fwhm_x = save_array[6]
                    self.fwhm_y = save_array[7]
                    self.finding_overlaps()
                    self.calculate_averages()
                    return list(self.plotting_fuda_data.keys())[0]
                except:
                    return False
        else:
            return False


    def fuda(self, data, peak_loc, name='prelim', ftol=1e-10):
        if len(data.shape)==2:
            data = numpy.expand_dims(data, axis=0)
        x = numpy.arange(len(data[0,0,:]))
        y = numpy.arange(len(data[0,:,0]))
        x,y = numpy.meshgrid(x,y)
        self.fuda_data = data
        self.resim_data = []
        intensities = []
        plotting_fuda_data_row = []
        plotting_resim_data_row = []
        fwhm_x_array = []
        fwhm_y_array = []
        amp_array = []

        y0_init, x0_init = peak_loc
        y0_init = int(y0_init)
        x0_init = int(x0_init)
        for number, dat in enumerate(data):
            initial_guess = []
            bounds_neg = []
            bounds_pos = []
            self.thresh=float(numpy.max(numpy.fabs(dat))*self.frac_thresh)

            if number == 0:
                initial_guess.append(self.thresh)
                bounds_neg.append(self.thresh-1e-10)
                bounds_pos.append(self.thresh+1e-10)

                amp_factor  = (numpy.reshape(PseudoVoigt((x,y), 0., 1., x0_init, y0_init, self.Gamma_x, self.Gamma_y, self.sigma_x, self.sigma_y, self.nu1, self.nu2), dat.shape)[y0_init, x0_init])
                initial_guess.append(dat[y0_init, x0_init]/amp_factor)
                bounds_neg.append(-numpy.inf)
                bounds_pos.append(numpy.inf)
                initial_guess.append(float(x0_init))
                bounds_neg.append(-numpy.inf)
                bounds_pos.append(numpy.inf)
                initial_guess.append(float(y0_init))
                bounds_neg.append(-numpy.inf)
                bounds_pos.append(numpy.inf)
                self.default_values_gauss(initial_guess, bounds_neg, bounds_pos)
                
                # print(numpy.max(numpy.reshape(PseudoVoigt((x,y), 0., 1., x0_init, y0_init, self.Gamma_x, self.Gamma_y, self.sigma_x, self.sigma_y, self.nu1, self.nu2), dat.shape)))
                # exit()
                initial_guess = list(initial_guess)
                print('maximum: ', numpy.max(dat*(numpy.fabs(dat)>self.thresh)).ravel())
                print(bounds_neg)
                print(bounds_pos)
                print(initial_guess)
                popt, pcov = opt.curve_fit(PseudoVoigt, (x,y), (dat*(numpy.fabs(dat)>self.thresh)).ravel(), p0 = initial_guess, verbose=1, maxfev=100000, ftol=ftol,gtol=1e-12,xtol=1e-14, bounds = (bounds_neg, bounds_pos))
                print('popt asdf', popt)
                self.save_fuda_values(popt)


            else:
                initial_guess.append(self.thresh)
                bounds_neg.append(self.thresh-1e-10)
                bounds_pos.append(self.thresh+1e-10)
                amp_factor  = (numpy.reshape(PseudoVoigt((x,y), 0., 1., x0_init, y0_init, self.Gamma_x, self.Gamma_y, self.sigma_x, self.sigma_y, self.nu1, self.nu2), dat.shape)[y0_init, x0_init])
                initial_guess.append(dat[y0_init, x0_init]/amp_factor)
                bounds_neg.append(-numpy.inf)
                bounds_pos.append(numpy.inf)
                initial_guess.append(x0_init)
                bounds_neg.append(x0_init-2.)
                bounds_pos.append(x0_init+2.)
                initial_guess.append(y0_init)
                bounds_neg.append(y0_init-2.)
                bounds_pos.append(y0_init+2.)
                initial_guess.append(self.Gamma_x)
                bounds_neg.append(self.Gamma_x-0.25*numpy.abs(self.Gamma_x))
                bounds_pos.append(self.Gamma_x+0.25*numpy.abs(self.Gamma_x))
                initial_guess.append(self.Gamma_y)
                bounds_neg.append(self.Gamma_y-0.25*numpy.abs(self.Gamma_y))
                bounds_pos.append(self.Gamma_y+0.25*numpy.abs(self.Gamma_y))
                initial_guess.append(self.sigma_x)
                bounds_neg.append(self.sigma_x-0.25*numpy.abs(self.sigma_x))
                bounds_pos.append(self.sigma_x+0.25*numpy.abs(self.sigma_x))
                initial_guess.append(self.sigma_y)
                bounds_neg.append(self.sigma_y-0.25*numpy.abs(self.sigma_y))
                bounds_pos.append(self.sigma_y+0.25*numpy.abs(self.sigma_y))
                initial_guess.append(self.nu1)
                bounds_neg.append(0)
                bounds_pos.append(1e-6)
                initial_guess.append(self.nu2)
                bounds_neg.append(0)
                bounds_pos.append(1e-6)
                initial_guess = list(initial_guess)
                popt, pcov = opt.curve_fit(PseudoVoigt, (x,y), (dat*(numpy.fabs(dat)>self.thresh)).ravel(), p0 = initial_guess, verbose=1, maxfev=100000, ftol=ftol,gtol=1e-12,xtol=1e-12, bounds = (bounds_neg, bounds_pos)) 
                self.save_fuda_values(popt)

            print(self.plotting_resim_data.keys())
            print(popt)
            # exit()
            resim_x = PV(1, numpy.linspace(0,data.shape[0], 200), data.shape[0]/2., popt[4], popt[6], popt[-2])[:101]
            resim_y = PV(1, numpy.linspace(0,data.shape[1], 200), data.shape[1]/2., popt[5], popt[7], popt[-1])[:101]
            frac_hwhm_x = 100-numpy.argmin(numpy.fabs(resim_x-(0.5*numpy.max(numpy.fabs(resim_x))))) # how far away are we from the centre of the PV function in simulated PV units?
            frac_hwhm_y = 100-numpy.argmin(numpy.fabs(resim_y-(0.5*numpy.max(numpy.fabs(resim_y))))) # how far away are we from the centre of the PV function in simulated PV units?
            hwhm_x=numpy.abs(data[0].shape[0]*frac_hwhm_x/200.) # convert between simulated PV units and actual data units
            hwhm_y=numpy.abs(data[0].shape[1]*frac_hwhm_y/200.) # convert between simulated PV units and actual data units
            fwhm_x = hwhm_y*2
            fwhm_y = hwhm_x*2

            fwhm_x = -self.uc1.ppm(fwhm_x)+ self.uc1.ppm(0)
            fwhm_y = -self.uc0.ppm(fwhm_y)+ self.uc0.ppm(0)

            fwhm_x_array.append(fwhm_x)
            fwhm_y_array.append(fwhm_y)

            amp = popt[1]
            x0 = popt[2]
            y0 = popt[3]

            amp_array.append(amp)



            # intensities.append(numpy.sum(dat))
            intensities.append(numpy.sum(PseudoVoigt((x,y),0,amp,x0,y0,self.Gamma_x, self.Gamma_y, self.sigma_x, self.sigma_y, self.nu1, self.nu2)))
            print('inten', intensities)
            self.resim_data.append(popt)
            plotting_fuda_data_row.append(dat)
            plotting_resim_data_row.append(popt)

        # try:
        #     self.plotting_resim_data[str(int(name))] = numpy.array(plotting_resim_data_row)
        #     self.plotting_fuda_data[str(int(name))] = numpy.array(plotting_fuda_data_row)

        # except:
        self.plotting_resim_data[str(name)] = numpy.array(plotting_resim_data_row)
        self.plotting_fuda_data[str(name)] = numpy.array(plotting_fuda_data_row)
        # self.overlap_data = 

        intensities = numpy.array(intensities)
        fwhm_x_array = numpy.array(fwhm_x_array)
        fwhm_y_array = numpy.array(fwhm_y_array)
        amp_array = numpy.array(amp_array)

        

        return intensities, fwhm_x_array, fwhm_y_array, amp_array

    def default_values(self, initial_guess, bounds_neg, bounds_pos):
        
        initial_guess.append(self.Gamma_x)
        bounds_neg.append(0)
        bounds_pos.append(200)
        initial_guess.append(self.Gamma_y)
        bounds_neg.append(0)
        bounds_pos.append(200)
        initial_guess.append(self.sigma_x)
        bounds_neg.append(0)
        bounds_pos.append(numpy.inf)
        initial_guess.append(self.sigma_y)
        bounds_neg.append(0)
        bounds_pos.append(numpy.inf)
        initial_guess.append(self.nu1)
        bounds_neg.append(0)
        bounds_pos.append(1)
        initial_guess.append(self.nu2)
        bounds_neg.append(0)
        bounds_pos.append(1)

    def default_values_gauss(self, initial_guess, bounds_neg, bounds_pos):
        
        initial_guess.append(self.Gamma_x)
        bounds_neg.append(0)
        bounds_pos.append(200)
        initial_guess.append(self.Gamma_y)
        bounds_neg.append(0)
        bounds_pos.append(200)
        initial_guess.append(self.sigma_x)
        bounds_neg.append(0)
        bounds_pos.append(numpy.inf)
        initial_guess.append(self.sigma_y)
        bounds_neg.append(0)
        bounds_pos.append(numpy.inf)
        initial_guess.append(0)
        bounds_neg.append(0.0)
        bounds_pos.append(1e-6)
        initial_guess.append(0)
        bounds_neg.append(0.)
        bounds_pos.append(1e-6)

        
