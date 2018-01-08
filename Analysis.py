"""
Analysis.py
Module to perform analysis on the data for the Oslo Model.
I Manco    20/02/17

Classes: - Measurements
         - Scaling
"""

import numpy as np
from matplotlib import pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy import stats
import Oslo_Model as om
reload(om)
import log_bin_CN_2016 as lb
reload(lb)

class Measurements:
    """Class that defines the properties of a given system."""
    
    def __init__(self, system_size, time, pile_height, avalanche_size, outflux):
    #def __init__(self, system_size, time, pile_height, avalanche_size):
        """Initialise a system with system size, # of iterations (= discrete time),
        pile height and avalanche size.
        """
        self.system_size = system_size
        self.time = time
        self.pile_height = pile_height
        self.avalanche_size = avalanche_size
        self.outflux = outflux
    
    def plot(self, time = "all"):
        """Plot pile height vs time for one system."""
        if time == "all":    
            plt.plot(self.time, self.pile_height)
        else:
            plt.plot(self.time[0:time], self.pile_height[0:time])
        plt.xlabel("Time")
        plt.ylabel("Pile Height (no. of grains)")
        plt.savefig("plot", dpi = 500)
    
    def moving_avg(self, w):
        moving_average = []
        time = []
        for t in range(0, self.time[-1], w):
            if t == 0:
                average = self.pile_height[0]
                moving_average.append(average)
                time.append(t)
            if t != 0:
                if t < w:
                    average = 1./t * sum(self.pile_height[0:2*t])
                    moving_average.append(average)
                    time.append(t)
                if t >= w:
                    average = 1./w * sum(self.pile_height[t-w:t+w])
                    moving_average.append(average)
                    time.append(t)
        return np.array(time), np.array(moving_average)
        
    def cross_over_time(self):
        """Return cross_over_time (depending on system size L)."""
        non_cross = []
        for i in range(1, self.time[-1], 10):
            if np.log(np.mean(self.pile_height[i+1000:i+2000])/np.mean(self.pile_height[i:i+1000])) >  0:
                non_cross.append(i)
            else:
                return i 

    def average_height(self):
        """Return average height after crossover time."""
        height = []
        for t in range(self.cross_over_time(), self.time[-1]):
            height.append(self.pile_height[t])
        return np.mean(height)
    
    def std_height(self):
        """Return the standard deviation of the height after the crossover time."""
        height_sqr = []
        for t in range(self.cross_over_time(), self.time[-1]):
            height_sqr.append(self.pile_height[t]**2)
        sqr_mean = np.mean(height_sqr)
        return np.sqrt(sqr_mean - self.average_height()**2)
    
    def height_prob(self):
        """Return the probability of the height after the crossover time."""
        heights = []
        for t in range(self.cross_over_time(), self.time[-1]):
            heights.append(self.pile_height[t])    
        x = np.bincount(heights)
        return x/float(self.time[-1])

    def avalanche_size_prob(self, N):
        """Task 3a. Only for system size L = 256"""
        avalanche = []
        if N == 1e6:
            for t in range(self.cross_over_time(), self.time[-1]):
                avalanche.append(self.avalanche_size[t])   
        else:
            for t in range(self.cross_over_time(), self.cross_over_time()+int(N)):
                avalanche.append(self.avalanche_size[t])      
        x = np.bincount(avalanche)
        return x/N

class Scaling:
    """Class to perform statistical analysis and scaling on a set of systems."""
    def __init__(self, N):
        N = int(N)
        system_a = om.results(8, N)
        system_b = om.results(16, N)
        system_c = om.results(32, N)
        #system_d = om.results(64, N)
        #system_e = om.results(128, N)
        #system_f = om.results(256, N)
        measurements = [system_a, system_b, system_c]
        #measurements = [system_a, system_b, system_c, system_d, system_e, system_f]
        self.measurements = measurements
        size = []
        cross_over_time = []
        for system in self.measurements:
            size.append(system.system_size)
            cross_over_time.append(system.cross_over_time())
        self.size = size
        self.cross_over_time = cross_over_time
        
    def mov_vs_time(self):
        """Plot moving average of height."""
        plt.figure()
        for system in self.measurements:
            plt.plot(np.log10(system.moving_avg(25)[0]), np.log10(system.moving_avg(25)[1]))
        plt.title("Log-log plot of moving average of height vs time")
        plt.xlim([0, 5])
        plt.xlabel("Log$_1$$_0$ (t)")
        plt.ylabel("Log$_1$$_0$ (h)")
        plt.show()
        
    def height_vs_time(self):
        """Task 2a. Plot height vs time for the various system sizes in the same plot."""
        plt.figure()
        for system in self.measurements:
            plt.loglog(system.time, system.pile_height, label = "L = " + str(system.system_size))
        plt.legend(loc = 2)
        plt.xlabel("$t$")
        plt.ylabel("$h$ ($t$; $L$)")
        plt.show()
    
    def time_vs_size(self):
        """Task 2a. Plot the cross-over-time vs system size and return fit."""
        def func(x, a, b):
            return a*(x**b)
        popt, pcov = curve_fit(func, self.size, self.cross_over_time)
        plt.figure()
        plt.plot(range(0, 260), popt[0]*(range(0, 260))**popt[1], '-', color = 'black')
        plt.plot(self.size, self.cross_over_time, 'x', color = 'r')
        plt.xlabel("$L$")
        plt.ylabel("$t$$_c$ ($L$)")
        plt.show()
        return popt, pcov
    
    def height_collapse(self):
        """Task 2b."""
        plt.figure()
        for system in self.measurements:
            x = system.moving_avg(100)[0]
            y = system.moving_avg(100)[1]
            #correction = 1.73*system.system_size**2 #instead of system.cross_over_time()
            plt.loglog(x/np.full(len(x), system.cross_over_time()), y/(np.full(len(y), system.system_size)), label = 'L = ' + str(system.system_size))
        plt.xlabel("$t$ / $t$$_c$")
        plt.ylabel("$\~h$ / $L$")
        plt.ylim(1.6*1e-1, 1e1)
        plt.legend(loc = "upper left")
        plt.show()

    def height_vs_size(self, a_0):
        """Task 2c. Plot average height in the steady state against the system size.
        Obtain scaling without correction.
        """
        size = []
        height = []
        for system in self.measurements:
            size.append(system.system_size)
            height.append(system.average_height())
        
        plt.plot(np.array(size), height, 'x', color = 'r')
        z = stats.linregress(size, height)
        sigma = z[-1]*np.sqrt(1./(sum(height)-np.mean(height))**2)
        plt.plot(range(260), np.full(260, z[0])*range(260) + z[1], color = 'black')
        plt.title("Average Height in Steady State")
        plt.xlabel("$L$")
        plt.ylabel("<$h$>")
        plt.figure()
        plt.plot(np.log10(size), np.array(height)/np.array(size), 'x', color = 'r', label = 'Without correction')
        plt.plot(np.log10(size), np.full(6, a_0), 'x', color = 'b', label = 'With correction')
        plt.xlabel("$L$")
        plt.ylabel("<$h$>/$L$")
        plt.legend(loc = 'lower right')
        plt.show()
        return z, sigma
    
    def correction(self):
        """Task 2c. Find correction to scaling."""
        size = []
        #height over size = average slope
        height_over_size = []
        for system in self.measurements:
            size.append(system.system_size)
            height_over_size.append(system.average_height()/system.system_size)
        plt.figure()
        plt.plot(size, height_over_size, '.')
        def func(x, a_0, a_1, omega):
            return a_0-a_0*a_1*(x**(-omega))
        popt, pcov = curve_fit(func, self.size, height_over_size)
        a_0 = popt[0]
        a_1 = popt[1]
        omega = popt[2]
        plt.plot(range(300), a_0 - a_0*a_1*range(300)**(-omega))
        return popt, pcov
        plt.title("Correction to scaling of average height with system size.")
        plt.xlabel("System Size L")
        plt.ylabel("<$h$>/$L$")
        plt.show()
      
    def find_a0(self):
        a = self.correction()[0][0]
        #omega = self.correction()[0][2]
        size = []
        height_over_size = []
        plt.figure()
        for system in self.measurements:
            size.append(system.system_size)
            height_over_size.append(system.average_height()/system.system_size)
        #for i in np.arange(a-.05, a+.01, .01):
        y = 1.73 - height_over_size
        x = size
        plt.loglog(x, y, '.', label = a)
        print x, y
        z = np.polyfit(np.array(x), np.array(y), 1)
        print z
        plt.xlabel('$L$')
        plt.ylabel('a$_0$ - <$h$>/$L$')
        plt.legend(loc = 'upper right')
        plt.show()
                  
    def std_vs_size(self):
        """Task 2c. Plot standard deviation of height/L vs system size."""
        size = []
        std = []
        for system in self.measurements:
            size.append(system.system_size)
            std.append(system.std_height())
        z =  np.polyfit(np.log10(size), np.log10(std), 1)
        plt.plot(np.log10(size), z[1]+ np.log10(size)*z[0], color = 'black')
        plt.plot(np.log10(size), np.log10(std), 'x', color = 'red')
        plt.figure()
        plt.loglog(size, std/size**z[0], 'x')
        #plt.plot((np.log10(size))**0.24, np.log10(std*size))
        plt.xlabel("$L$")
        plt.ylabel("$\sigma$$_h$/$L$")
        plt.legend()
        plt.show()
        return z
    
    def prob_vs_height(self):
        """Task 2d. Plot Probability of height vs height for various system sizes on the same plot."""
        plt.figure()
        for system in self.measurements:
            x = np.arange(0, np.max(system.pile_height)+1)
            y = system.height_prob()
            plt.plot(x, y, label =  "L = " + str(system.system_size))
        plt.xlabel("$h$")
        plt.ylabel("P$_N$ ($h$; $L$)")
        plt.legend()
        plt.show()

    def prob_collapse(self):
        """Task 2d. Data collapse for the height probability."""
        plt.figure()
        for system in self.measurements:
            interval = np.arange(0, np.max(system.pile_height)+1)
            prob = system.height_prob()
            x = (interval - system.average_height())/system.std_height()
            y = (prob)*(system.std_height())
            plt.plot(x, y, label = '$L$ = '+ str(system.system_size))
        plt.xlabel("($h$ - <$h$>)/$\\sigma$")
        plt.ylabel("$\sigma$$_h$ P$_N$ ($h$; $L$)")
        plt.xlim(-6, 6)
        plt.legend()
        plt.show()

    def prob_vs_avasize(self, a = 1.4):
        """Task 3a. Plot raw probability of avalanche size vs avalanche size for L = 256"""
        system = self.measurements[-1]
        for N in [1e4, 1e5, 1e6]:
            x = np.arange(0, len(system.avalanche_size_prob(N)))
            y = system.avalanche_size_prob(N)
            if N == 1e4:
                lab = '$N$ = 10$^4$'
            if N == 1e5:
                lab = '$N$ = 10$^5$'
            if N ==1e6:
                lab = '$N$ = 10$^6$'
            plt.loglog(x, y, '.', label = lab, markersize = 3.)
        x = self.log_bin_data(system, a)[0]
        y = self.log_bin_data(system, a)[1]
        plt.loglog(x, y, '-', label = "$\~P$$_N$($s$; $L$)", color = 'black', linewidth = 1.)
        plt.title('$L$ = 256')
        plt.xlabel("$s$", fontsize = 16)
        plt.ylabel("P$_N$ ($s$; $L$)", fontsize = 14)
        plt.legend()
        plt.show()
  
    def log_bin_data(self, system, a = 1.4):
        os.chdir('/Users/Ilaria/Documents/Imperial/Third Year/C&N/Data/with_avalanche_size/' + str(int(1e6)) +'/'+str(system.system_size))
        f = open('Avalanche Size.txt', 'r')
        x = []
        for line in f:
            x.append(float(line))
        os.chdir('/Users/Ilaria/Documents/Imperial/Third Year/C&N/Code')
        x = np.array([int(z) for z in x])
        vals, counts = lb.frequency(x)
        b, c = lb.log_bin(x, 1., 1.5, a, 'integer', debug_mode=False, drop_zeros = False)
        """vals, counts = lb.lin_bin(x, int(max(x)))
        b, c = lb.log_bin(x, 1., 1., a, debug_mode=False)"""
        #return the binned data for each system size
        return np.array(b), np.array(c)
        
    def find_tau(self, a = 1.4):
        """Extrapolate avalanche-size exponent from slope of the PDF."""
        x = []
        system = self.measurements[-1]
        index = 0
        for i in self.log_bin_data(system, a)[0]:
            if i > 1e1:
                if i < 1e4:
                    x.append(i)
            else:
                index += 1
        c = self.log_bin_data(system, a)[1]
        #z =  np.polyfit(np.log10(x), np.log10(c[:len(x)]), 1)
        z = stats.linregress(np.log(x), np.log(c[index:len(x)+index]))
        sigma = z[-1]*np.sqrt(1./(sum(np.log(x)-np.log(np.mean(x))))**2)
        return z, sigma
        
    def plot_bin_data(self, a = 1.4):
        for system in self.measurements:
            x = self.log_bin_data(system, a)[0]
            y = self.log_bin_data(system, a)[1]
            plt.loglog(x, y, '-', label = '$L$ = '+str(system.system_size))
        plt.loglog(x, x**-1.55, '--', color = 'black', label = '$s$$^-$$^1$$^.$$^5$$^5$')
        plt.xlabel('s')
        plt.ylabel('$\~P$$_N$ (s)')
        plt.legend()
        plt.show()

    def collapse_ava(self, D = 0):
        """Plot collapse of the avalanche size probability."""
        plt.figure()
        tau = self.find_tau()[0][0]
        for system in self.measurements:
            x = (self.log_bin_data(system)[0])/(system.system_size**D)
            #x = self.log_bin_data(system)[0]/self.find_D()[0] #only y axis collapsed
            y = self.log_bin_data(system)[1]*((self.log_bin_data(system)[0])**(-tau))
            plt.loglog(x, y, label = '$L$ = '+str(system.system_size))
            """for i, j in enumerate(x):
                if j == max(x):
                    a = i 
            s_c.append(max(x))
            plt.loglog(max(x), y[a], 'x', color = 'black', markersize = 4.)"""
        plt.legend(loc = 'lower left')
        plt.xlabel("s/L$^D$")
        #plt.xlabel('$s$/$s$$_c$')
        plt.ylabel('$\~P$$_N$ (s)  s$^$\tau$')
        plt.show()
    
    def find_D(self):
        s_c = []
        for system in self.measurements:
            x = self.log_bin_data(system)[0] #only y axis collapsed
            s_c.append(max(x))
        #plt.loglog(np.log([8, 16, 32, 64, 128, 256]), np.log(s_c), 'x', color = 'red')
        D = stats.linregress(np.log([8, 16, 32, 64, 128, 256]), np.log(s_c))
        #plt.loglog(np.log([8, 16, 32, 64, 128, 256]), np.log([8, 16, 32, 64, 128, 256])*D[0], color = 'black')
        #plt.ylabel('$s$$_c$')
        #plt.xlabel('$L$')
        #plt.show()
        return s_c, D
    
    def moment(self, system, k):
        """Return kth moment for a given system."""
        s_k = []
        T = system.time[-1] - system.cross_over_time()
        for t in range(system.cross_over_time(), system.time[-1]):
            s_k.append(system.avalanche_size[t]**k)   
        moment = ((1./T))*sum(s_k)
        return moment
    
    def moment_scaling(self, first):
        grad = []
        plt.figure()
        for k in range(1, 6):
            size = []
            mom = []
            for system in self.measurements:
                size.append(system.system_size)
                mom.append(self.moment(system, k))
            mom = np.array(mom)
            size = np.array(size)
            plt.plot(np.log(size), np.log(mom), 'x', label = '$k$ = ' + str(k))
            slope, interc, rvalue, pvalue, stderr = stats.linregress(np.log(size)[first:], np.log(mom)[first:])
            grad.append(slope)
            plt.plot(np.log(size), slope*np.log(size) + interc)
            print slope
        plt.xlabel('System Size $L$')
        plt.ylabel('<s$^k$>')
        plt.legend(loc = 'upper left')
        plt.show()
        plt.figure()
        plt.plot(range(1, 6), grad, 'x', color = 'black')
        plt.xlabel("$k$")
        plt.ylabel("$D$ ($k$ - $\\tau$ + 1)")
        plt.xlim(0, 6)
        slope1, interc1, rvalue1, pvalue1, stderr1 = stats.linregress(range(1, 6), grad)
        plt.plot(range(1, 6), slope1*(np.array(range(1, 6))) + interc1, color = 'black')
        error1 = slope1*np.sqrt(1./(sum(np.log(grad)-np.log(np.mean(grad))))**2)
        return slope1, interc1, error1
    
    def outflux_prob(self):
        """Bonus Task 1."""
        for system in self.measurements[0:2]:
            #exclude d = 0
            outflux = []
            for d in system.outflux:
                if d != 0:
                    outflux.append(d)
            plt.hist(sorted(outflux), bins = len(outflux)*system.system_size*200, label = "$L$ = " + str(system.system_size))
        plt.xlabel('$d$')
        plt.legend()
        plt.show()