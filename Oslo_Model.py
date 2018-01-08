"""
Oslo_Model.py
Module to implement the algorithm for the Oslo Model.
I Manco    20/02/17

Classes: - System
         - Relaxation

Functions: - results(system_size, N, directory)

To save data in a directory, type path of desired folder into argument.
"""

import numpy as np
import random as rn
import Analysis as an
import os
#reload(an)

class System:
    """Initialisation of the system."""
    def __init__(self, L):
        self.L = L        
        self.sites = np.array(range(1, L + 1))
        self.h = np.full(len(self.sites), 1, int)
        threshold = []
        for i in range(0, len(self.sites)):
            ran = rn.randrange(1, 3)
            ran = np.random.choice([1, 2], p=[0.5, 0.5])
            threshold.append(ran)
        self.threshold = threshold
        
    def slope(self):
        h = self.h
        z = []
        for i in range(0, len(h)):
            if i == len(h) - 1:
                z.append(h[i])
                #z.append(1)
            else:
                slope = h[i] - h[i + 1] 
                z.append(slope)
        return z
    
    def total_height(self):
        return self.h[0]
                

class Relaxation:
    """Drive and relaxation of the system."""
    def __init__(self, system):
        self.system = system
    
    def drive(self):
        system = self.system
        system.h[0] = system.h[0] + 1
    
    def relax_site(self, i):
        system = self.system
        if i < len(system.sites) - 1:
            system.h[i] = system.h[i] - 1
            system.h[i + 1] = system.h[i + 1] + 1
            system.threshold[i] = np.random.choice([1, 2], p=[0.5, 0.5])
        if i == len(system.sites) -1 :
            system.h[i] = system.h[i] - 1
            system.threshold[i] = np.random.choice([1, 2], p=[0.5, 0.5])
            
    def check_relax(self):
        system = self.system
        for i in range(len(system.sites)):
            if system.slope()[i] > system.threshold[i]:
                return False

    def relax_system(self):
        system = self.system
        self.drive()
        avalanche_size = 0
        dropsize = 0
        while self.check_relax() ==  False:
            for i in range(len(system.sites)):
                if system.slope()[i] > system.threshold[i]:
                    self.relax_site(i) 
                    avalanche_size += 1
                    if i == len(system.sites) -1 :
                        dropsize += 1

        system.h
        return avalanche_size, dropsize
    
    def iterate(self, iterations):
        height = []
        time = []
        system = self.system
        avalanche_size = []
        dropsize = []
        for i in range(1, iterations):
            avalanche_size.append(self.relax_system()[0])
            height.append(system.total_height())
            time.append(i)
            if dropsize != 0:
                dropsize.append(self.relax_system()[1])
        return time, height, avalanche_size, dropsize
        

def results(system_size, N, directory):
    """Function that reproduces the ricepile experiment for different system sizes L and stores 
    the results in a text file if the experiment has not been carried out before. If data for the chosen 
    system size L already exists, the function accesses the files with the results.
    Also creates an instance of the the class Measurements to analyse the results."""
    
    system = System(system_size)
    relaxation = Relaxation(system)
    newdir = 'path'+str(int(N))+'/'+str(system_size)
    #check if directory with results for system size = L already exists (no need for experiment)
    if not os.path.exists(newdir):
        experiment = relaxation.iterate(N)
        times = experiment[0]
        heights = experiment[1]
        avalanche_sizes = experiment[2]
        dropsize = experiment[3]

        os.makedirs(newdir)
        os.chdir(str(directory)+str(system_size))
        f = open('Pile Height', 'w')
        f1 = open('Time', 'w')
        f2 = open('Avalanche Size', 'w')
        f3 = open('Drop Size', 'w')
        for height in heights:
            f.write(str(height) + '\n')
        for time in times:
            f1.write(str(time) + '\n')
        for ava_size in avalanche_sizes:
            f2.write(str(ava_size) + '\n')
        for drop in dropsize:
            f3.write(str(drop) + '\n')  
    os.chdir(str(directory)+str(int(N))+'/'+str(system_size))     
    f = open('Pile Height.txt', 'r')
    f1 = open('Time.txt', 'r')
    f2 = open('Avalanche Size.txt', 'r')
    f3 = open('Drop Size', 'r')
    time = []
    height = []
    ava_size = []
    outflux = []
    for line in f1:
        time.append(int(line))
    for line in f:
        height.append(int(line))
    for line in f2:
        ava_size.append(int(line))
    for line in f3:
        outflux.append(int(line))
    #results = an.Measurements(system_size, time, height, ava_size)
    results = an.Measurements(system_size, time, height, ava_size, outflux)
    os.chdir(str(directory)
    return results      

def outflux(system_size, N, directory):
    system = System(system_size)
    relaxation = Relaxation(system)
    newdir = str(directory)+str(int(N))+'/'+str(system_size)+'dropsize'
    #check if directory with results for system size = L already exists (no need for experiment)
    if not os.path.exists(newdir):
        experiment = relaxation.iterate(N)
        dropsize = experiment[3]
        os.chdir(str(directory)+str(int(N))+'/'+str(system_size))
        f = open('Drop Size', 'w')
        for drop in dropsize:
            f.write(str(drop) + '\n')  
    outflux = []
    f = open('Drop Size', 'r')
    for line in f:
        outflux.append(int(line))
    return outflux
