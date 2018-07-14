# -*- coding: utf-8 -*-
"""
Created on Tue May 15 10:46:14 2018

@author: c1341133

Calibration

"""
import numpy as np
import xml.etree.ElementTree as ET
import pyfits

def loadFits(fileName):
    hdul = pyfits.open(fileName)
    prihdr = hdul[0].header
    
    
    NDR0_time = hdul[-2].data['Time 1.0 s'][0]
    NDRfinal_time = hdul[-2].data['Time 1.0 s'][prihdr['MACCUM']-1]
    CDS_time =  NDRfinal_time - NDR0_time

    info = {'NEXP': prihdr['NEXP'], 'MACCUM': prihdr['MACCUM'], 'TEXP': prihdr['TEXP'],
            'WL': hdul[-4].data['Wavelength 1.0 um'], 'Z': hdul[-2].data['z'],
            'CDS_time': CDS_time, 'T': hdul[-2].data['Time 1.0 s'] }
        
    data = np.zeros((hdul[1].data.shape[0], hdul[1].data.shape[1],info['NEXP']*info['MACCUM']))
    for i in range(data.shape[2]):
        data[...,i] = hdul[1+i].data
    
    return data, info
    
    

def doCDS(data, info):
    
    cds_data = np.zeros(( data.shape[0], data.shape[1] ,info['NEXP']))

#    for i in range (info['NEXP']):
    for i in range (data.shape[2]/2):


        idx_1 =  i*info['MACCUM']   
        idx_2 =  idx_1 + (info['MACCUM']  - 1)
        cds_data[...,i] = data[...,idx_2] - data[...,idx_1]
    

    return cds_data
            


def subtractDark(data, info, ICF, ch):

    root = ET.parse(ICF).getroot()
    for child in root.findall('channel'):
        if child.get('name') == ch:
            dc = np.float(child.find('detector_pixel').find('Idc').get('val'))
       
    data = data - dc*info['CDS_time']
    
    return data
    
   
def flatField(data, ICF, ch, QE_rms_file):

    QE_grid = np.load(QE_rms_file)[:data.shape[0],:data.shape[1]]
    
    data =  np.rollaxis(data,2,0)
    data = data/QE_grid
    data =  np.rollaxis(data,0,3) 
    
    return data
            

def backSub(data):
    
    for i in range (data.shape[2]):
   
       background_1 = data[...,i][5:10,10:-10]
       background_2 = data[...,i][-10:-5,10:-10]   
       background = np.mean(np.vstack ( (background_1, background_2) ) )
       data[...,i] -= background
    
    return data
        
 