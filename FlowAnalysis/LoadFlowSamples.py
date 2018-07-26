# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 15:22:47 2016
Tools for flow cytometry analysis. Bases on FlowCytometryTools(https://pypi.python.org/pypi/FlowCytometryTools/0.4.0)
Useful for dealing with multiple samples at once.

Example:

Loading multiple samples by pattern:
import LoadFlowSamples_v1 as LFS
pattern = 'Tube'
datadir = '/Users/Alonyan/Documents/Python/Comps'
LFS.FlowData(datadir, pattern).samples

Gating: 
Samples = LFS.gate(Samples, gate1)
gate1 is a gate object f.i.
from FlowCytometryTools import PolyGate
gate1 = PolyGate([(3.393e+04, 6.062e+04), (4.619e+04, 1.696e+04), (1.306e+05, 2.771e+04), (2.311e+05, 1.050e+05), (2.538e+05, 1.829e+05), (2.090e+05, 2.393e+05), (1.656e+05, 2.440e+05), (9.007e+04, 2.091e+05), (5.280e+04, 1.446e+05), (5.280e+04, 1.446e+05), (5.280e+04, 1.446e+05)], ('FSC-A', 'SSC-A'), region='in', name='gate1')

Applying transformations (['hlog' | 'tlog' | 'glog' | callable])
Trans_Samples = LFS.transform(Samples, 'hlog', Channels, 500)

Applying compensations:
NewSamps = LFS.compansate(CompensationSamples, CompMat, Channels)

Calculate summary stats
LFS.means(Samples, ['Alexa Fluor 488-A'])

@author: Alonyan
"""

from FlowCytometryTools import FCMeasurement
import os
import numpy as np
import pandas as pd


class FlowData:
    def __init__(self, datadir, pattern):
        self.datadir = datadir
        self.pattern = pattern
        self.samples = []
        self.samples = self._load()
                
        blankInd = []
        if 'FSC-A' in self.samples[0].channel_names:
            blankInd.append(self.samples[0].channel_names.index('FSC-A'))
        if 'FSC-H' in self.samples[0].channel_names:
            blankInd.append(self.samples[0].channel_names.index('FSC-H'))
        if 'FSC-W' in self.samples[0].channel_names:
            blankInd.append(self.samples[0].channel_names.index('FSC-W'))
        if 'SSC-A' in self.samples[0].channel_names:
            blankInd.append(self.samples[0].channel_names.index('SSC-A'))
        if 'SSC-H' in self.samples[0].channel_names:
            blankInd.append(self.samples[0].channel_names.index('SSC-H'))
        if 'SSC-W' in self.samples[0].channel_names:
            blankInd.append(self.samples[0].channel_names.index('SSC-W'))
        if 'Time' in self.samples[0].channel_names:
            blankInd.append(self.samples[0].channel_names.index('Time'))
            
        ChanInd = list(set(np.arange(0,len(self.samples[0].channel_names)))-set(blankInd))
        ChanList = list(np.array(self.samples[0].channel_names)[ChanInd])
        self.fluorescent_channel_names = list(ChanList)
        self.channel_names = list(self.samples[0].channel_names)
        self.samples = _asinhtform(self.samples, self.fluorescent_channel_names, 150)
        self.IDlist = [self.samples[0].ID[0:-4]]
        for i in np.arange(1,len(self.samples)):
            self.IDlist = self.IDlist + [self.samples[i].ID[0:-4]]

        

        
    def _load(self):
        for file in os.listdir(self.datadir):
            if file.endswith(".fcs") and self.pattern in file:
                self.samples.append(FCMeasurement(ID=file, datafile=self.datadir+"/"+ file))
        return np.array(self.samples)



    
    def copy(self):
        import copy
        return copy.deepcopy(self)
    def LinearRange(self):
        return 150
        
    def means(self):
        return _means(self.samples, self.fluorescent_channel_names)
    def aSinhmeans(self):#means in arcsinh space. Use for compensations
        return _aSinhmeans(self.samples, self.fluorescent_channel_names)
    def medians(self):
        return _medians(self.samples, self.fluorescent_channel_names)
    def stds(self):
        return _stds(self.samples, self.fluorescent_channel_names)
    def counts(self):
        return _counts(self.samples, [self.channel_names[0]])
    def geoMeans(self):
        return _geoMeans(self.samples, self.fluorescent_channel_names)
    
    def compensate(self, CompMat):
        Comped = self.copy()
        Comped.samples = _compansate(Comped.samples, CompMat,Comped.fluorescent_channel_names)
        return Comped
    
    def gate(self, gate1):
        gated = self.copy()
        gated.samples = _gate(gated.samples, gate1)
        return gated
    
    def transform(self, tform, Channels, LinRange):
        transformed = self.copy()
        transformed.samples = _transform(transformed.samples, tform, Channels, LinRange)
        return transformed
    
#    def asinhtform(self, Channel, LinRange):
#        transformed = self.copy()
#        transformed.samples = _asinhtform(transformed.samples, Channel, LinRange)
#        return transformed




      
def _gate(samples, gate1): 
    Gatedsamples = samples
    for i in np.arange(0,len(samples)):
        Gatedsamples[i] = samples[i].gate(gate1)
    return Gatedsamples
    
def _transform(samples, tform, Channels, LinRange):
    Transsamples = samples
    for i in np.arange(0,len(samples)):
        Transsamples[i] = samples[i].transform(tform, channels = Channels, b=LinRange)
    return Transsamples
    
def _asinhtform(samples, Channel, LinRange):
    Transsamples = samples
    i=0
    for sample in samples: 
        tmpsample = sample.copy()
        new_data = tmpsample.data
        new_data[Channel] = np.arcsinh(sample.data[Channel]/LinRange)
        tmpsample.data = new_data
        Transsamples[i]=tmpsample
        i+=1
    return Transsamples

def _compansate(samples, CompMat, Channels):
    CompSamples = samples.copy()
    i=0
    for sample in samples: 
        tmpsample = sample.copy()
        new_data = tmpsample.data
        new_data[Channels] = CompMat.dot(sample.data[Channels].T).T
        tmpsample.data = new_data
        CompSamples[i]=tmpsample
        i+=1
    return CompSamples


def _aSinhmeans(samples,ChanList):
    if type(samples) is np.ndarray:       
        aSinhmeans = samples[0][ChanList].mean()      
        IDlist = [samples[0].ID[0:-4]]
        for i in np.arange(1,len(samples)):
            aSinhmeans = pd.concat([aSinhmeans, samples[i][ChanList].mean()], axis=1, ignore_index=True)
            IDlist = IDlist + [samples[i].ID[0:-4]]
    else:
        aSinhmeans = samples[ChanList].mean()      
        IDlist = [samples.ID[0:-4]]
    aSinhmeans.columns = IDlist
    return aSinhmeans.T
#    
#    
    
def _means(samples,ChanList):
    if type(samples) is np.ndarray:
        means = np.mean(150*np.sinh(samples[0][ChanList]))
        IDlist = [samples[0].ID[0:-4]]
        for i in np.arange(1,len(samples)):
            means = pd.concat([means, np.mean(150*np.sinh(samples[i][ChanList]))], axis=1, ignore_index=True)
            IDlist = IDlist + [samples[i].ID[0:-4]]
    else:
        means = np.mean(150*np.sinh(samples[ChanList]))
        IDlist = [samples.ID[0:-4]]
    means.columns = IDlist
    return means.T
        
def _medians(samples,ChanList):
    if type(samples) is np.ndarray:  
        medians = 150*np.sinh(samples[0][ChanList].median())      
        IDlist = [samples[0].ID[0:-4]]
        for i in np.arange(1,len(samples)):
            medians = pd.concat([medians, 150*np.sinh(samples[i][ChanList].median())], axis=1, ignore_index=True)
            IDlist = IDlist + [samples[i].ID[0:-4]]
    else:
        medians = 150*np.sinh(samples[ChanList].median())      
        IDlist = [samples.ID[0:-4]]
    medians.columns = IDlist
    return medians.T
        
#def _stds(samples,ChanList):
#    if type(samples) is np.ndarray:
#        stds = samples[0][ChanList].std()      
#        IDlist = [samples[0].ID[0:-4]]
#        for i in np.arange(1,len(samples)):
#            stds = pd.concat([stds, samples[i][ChanList].std()], axis=1, ignore_index=True)
#            IDlist = IDlist + [samples[i].ID[0:-4]]
#    else:
#        stds = samples[ChanList].std()      
#        IDlist = [samples.ID[0:-4]]
#    stds.columns = IDlist
#    return stds.T        

def _stds(samples,ChanList):
    if type(samples) is np.ndarray:
        stds = np.std(150*np.sinh(samples[0][ChanList]))
        IDlist = [samples[0].ID[0:-4]]
        for i in np.arange(1,len(samples)):
            stds = pd.concat([stds, np.std(150*np.sinh(samples[i][ChanList]))], axis=1, ignore_index=True)
            IDlist = IDlist + [samples[i].ID[0:-4]]
    else:
        stds = np.std(150*np.sinh(samples[ChanList]))
        IDlist = [samples.ID[0:-4]]
    stds.columns = IDlist
    return stds.T       
        
        
def _counts(samples,ChanList):
    if type(samples) is np.ndarray:
        counts = samples[0][ChanList].count()      
        IDlist = [samples[0].ID[0:-4]]
        for i in np.arange(1,len(samples)):
            counts = pd.concat([counts, samples[i][ChanList].count()], axis=1, ignore_index=True)
            IDlist = IDlist + [samples[i].ID[0:-4]]
    else:
        counts = samples[ChanList].count()      
        IDlist = [samples.ID[0:-4]]
    counts.columns = IDlist
    counts.index = ['# Cells']
    return counts.T
        
        
 
def _geoMeans(samples,ChanList):
    geoMeans = 150*np.sinh(_aSinhmeans(samples,ChanList))
    return geoMeans