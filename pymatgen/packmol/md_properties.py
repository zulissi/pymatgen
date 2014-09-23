# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 01:16:09 2013

@author: 
"""

import lammpsio
import sys
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
from multiprocessing import Pool
from pymatgen.packmol.lammpsio import LammpsLog


def autocorrelate (a):
    b=np.concatenate((a,np.zeros(len(a))),axis=1)
    c= np.fft.ifft(np.fft.fft(b)*np.conjugate(np.fft.fft(b))).real
    d=c[:len(c)/2]
    d=d/(np.array(range(len(a)))+1)[::-1]
    return d[:wline]

if __name__=='__main__':
    
    """
    All the properties are evaluated based on the *properties input. Right now, pymatgen only supports viscosity (argument: viscosity) evaluation.
                Example: 'md_properties log.lammps viscosity' will return the viscosity of the system.
    Modified on June 11, 2014 by YZ: read window size and cutoff from command line
                                     window size: # of lines of the window
                                     cutoff: initial lines ignored for the calculation
                Example: md_properties log.lammps 4 0
    """
    
    logfilename=sys.argv[1]
    wline=int(sys.argv[2])
    cutoff=float(sys.argv[3])
    properties=sys.argv[4:]
    l=LammpsLog.from_file(logfilename)
    #wline=l.wline
    #print wline
    print 'Done reading the log file. Starting Calculations...'
    NCORES=8
    p=Pool(NCORES)
    
  

#    if 'viscosity' in l.properties:
    a1=l.llog['pxy']
    a2=l.llog['pxz']
    a3=l.llog['pyz']
    a4=l.llog['pxx']-l.llog['pyy']
    a5=l.llog['pyy']-l.llog['pzz']
    a6=l.llog['pxx']-l.llog['pzz']
    array_array=[a1,a2,a3,a4,a5,a6]
    pv=p.map(autocorrelate,array_array)
    pcorr = (pv[0]+pv[1]+pv[2])/6+(pv[3]+pv[4]+pv[5])/24
    
    
    print l.llog['timestep']
    temp=np.mean(l.llog['temp'][cutoff:])
    print temp
    print cutoff
    visco = (scipy.integrate.cumtrapz(pcorr,l.llog['step'][:len(pcorr)]))*l.llog['timestep']*10**-15*1000*101325.**2*l.llog['vol'][-1]*10**-30/(1.38*10**-23*temp)
    plt.plot(np.array(l.llog['step'][:len(pcorr)-1])*l.llog['timestep'],visco)
    plt.xlabel('Time (femtoseconds)')
    plt.ylabel('Viscosity (cp)')
    plt.savefig('viscosity_parallel.png')
    
    output=open('viscosity_parallel.txt','w')
    output.write('#Time (fs), Average Pressure Correlation (atm^2), Viscosity (cp)\n')
    for line in zip(np.array(l.llog['step'][:len(pcorr)-1])*l.llog['timestep']-cutoff,pcorr,visco):
        output.write(' '.join(str(x) for x in line)+'\n')
    output.close()
    print 'Viscosity Calculation Comlete!'
