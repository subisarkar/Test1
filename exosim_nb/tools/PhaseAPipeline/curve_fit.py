# -*- coding: utf-8 -*-
"""
Created on Tue May 29 07:28:53 2018

@author: user1
"""

import numpy as np
import pytransit
from scipy.optimize import minimize
import matplotlib.pyplot as plt


modelMandelAgol =  pytransit.MandelAgol(eclipse=False)

def lc_model(z, p, gamma, multiaccum, NDR0_time, NDRfinal_time):    
    LC = modelMandelAgol(z , np.sqrt(p), gamma)    
    LC0 = LC[0: len(LC) :multiaccum] * NDR0_time
    LC1 = LC[multiaccum-1: len(LC) :multiaccum] * NDRfinal_time    
    LC = LC1-LC0
    LC = LC/ LC[0]
    return LC
    
def chi_sq ((p, F, g1, g2), lc, err, z, multiaccum, NDR0_time, NDRfinal_time):
    model = lc_model(z, p, [g1,g2],  multiaccum, NDR0_time, NDRfinal_time)  * F    
    return np.sum(((model-lc))**2)

def light_curve_fit(ex, z, multiaccum, NDR0_time, NDRfinal_time): 
    err = np.std(np.hstack((ex[0:np.int(0.2*len(ex))], ex[np.int(0.8*len(ex)):])))
    oot_est = np.mean(np.hstack((ex[0:np.int(0.2*len(ex))], ex[np.int(0.8*len(ex)):])))
    it_est = np.mean(ex[len(ex)/2 -len(ex)/8: len(ex)/2 +len(ex)/8])
    cr_est = (oot_est - it_est ) / oot_est
    if err == 0:
        err = 1e-8
              
    fit_init = [cr_est, oot_est, 0.0, 0.0]
    fit  = minimize(chi_sq,fit_init, args=(ex, err, z, multiaccum, NDR0_time, NDRfinal_time), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=None, callback=None, options=None) 
         
    return fit['x']

def find_transit_depth(data, info):
    
#    data =  np.load('/Users/user1/ExoSimOutput_light_curve/SvRuns/GJ 1214 b_NIR Spec_sim0000.npy')
#    z = np.load('/Users/user1/ExoSimOutput_light_curve/SvRuns/GJ 1214 b_NIR Spec_sim0000-z.npy')
#    t = np.load('/Users/user1/ExoSimOutput_light_curve/SvRuns/GJ 1214 b_NIR Spec_sim0000-t.npy')
  
    z = info['Z']
    multiaccum = info['MULTIACCUM']
    t = info['T']
    
    NDR0_time = t[0]
    NDRfinal_time = t[2-1]
    multiaccum = 2   
    
    transit_depths = np.zeros((data.shape[1]))
    err_depths = np.zeros((data.shape[1]))

    for i in range (data.shape[1]) :
        
        lc = data[:,i]
        fit = light_curve_fit(lc, z, multiaccum, NDR0_time, NDRfinal_time)

        lc_m = lc_model(z, fit[0], [fit[2], fit[3]], multiaccum, NDR0_time, NDRfinal_time)
        
        residual = (lc/lc.max()) - lc_m
        
#        plt.figure('lc')
#        plt.plot(lc_m)
#        plt.plot(lc/lc.max(), 'rx')
#        plt.figure('res')
#        plt.plot(residual)
             
        transit_depths[i] = fit[0]
        err_depths[i] = np.sqrt(2)* np.std(residual)/np.sqrt(len(lc/2))
        
    return transit_depths, err_depths


