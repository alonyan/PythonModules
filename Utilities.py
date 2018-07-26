# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 15:46:30 2016

@author: Alonyan
"""
import numpy as np
import scipy.io as sio

def num2str(num, precision): 
    return "%0.*f" % (precision, num)
    
def colorcode(datax, datay):
    from scipy import interpolate
    import numpy as np
    H, xedges, yedges = np.histogram2d(datax,datay, bins=30)
    xedges = (xedges[:-1]+xedges[1:])/2
    yedges = (yedges[:-1]+yedges[1:])/2
    f = interpolate.RectBivariateSpline(xedges,yedges , H)
    
    z = np.array([])
    for i in datax.index:
        z = np.append(z,f(datax[i],datay[i]))
    #z=(z-min(z))/(max(z)-min(z))   
    z[z<0] = 0
    idx = z.argsort()
    return z, idx
    
class kmeans:        
    def __init__(self, X, K):
        # Initialize to K random centers
        oldmu = X.sample(K).values#np.random.sample(X, K)
        mu = X.sample(K).values#np.random.sample(X, K)
        while not _has_converged(mu, oldmu):
            oldmu = mu
            # Assign all points in X to clusters
            clusters = _cluster_points(X, mu)
            # Reevaluate centers
            mu = _reevaluate_centers(oldmu, clusters)
        self.mu = mu
        self.clusters = clusters
        #return(mu, clusters)
        
def _cluster_points(X, mu):
    clusters  = {}
    for x in X:
        bestmukey = min([(i[0], np.linalg.norm(x-mu[i[0]])) \
                    for i in enumerate(mu)], key=lambda t:t[1])[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    return clusters
 
def _reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis = 0))
    return newmu
 
def _has_converged(mu, oldmu):
    return (set(mu) == set(oldmu))
    
    
    

#Utils for opening MAT files
def print_mat_nested(d, indent=0, nkeys=0):
    """Pretty print nested structures from .mat files   
    Inspired by: `StackOverflow <http://stackoverflow.com/questions/3229419/pretty-printing-nested-dictionaries-in-python>`_
    """
    
    # Subset dictionary to limit keys to print.  Only works on first level
    if nkeys>0:
        d = {k: d[k] for k in d.keys()[:nkeys]}  # Dictionary comprehension: limit to first nkeys keys.

    if isinstance(d, dict):
        for key, value in d.iteritems():         # iteritems loops through key, value pairs
          print('\t' * indent + 'Key: ' + str(key))
          print_mat_nested(value, indent+1)

    if isinstance(d,np.ndarray) and d.dtype.names is not None:  # Note: and short-circuits by default
        for n in d.dtype.names:    # This means it's a struct, it's bit of a kludge test.
            print('\t' * indent + 'Field: ' + str(n))
            print_mat_nested(d[n], indent+1)


def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    
    from: `StackOverflow <http://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries>`_
    '''
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    dict1 = {}
    for key in dict:
        if isinstance(dict[key],np.ndarray):
            i=1
            for inst in dict[key]:
                if isinstance(inst, sio.matlab.mio5_params.mat_struct):
                    dict1[key+'_'+str(i)] = _todict(inst) 
                    i+=1
        elif isinstance(dict[key], sio.matlab.mio5_params.mat_struct):
            dict1[key] = _todict(dict[key])
            
    return dict1       

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, sio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        elif isinstance(elem,np.ndarray):
            dict[strg] = _tolist(elem)
        else:
            dict[strg] = elem
    return dict

def _tolist(ndarray):
    '''
    A recursive function which constructs lists from cellarrays 
    (which are loaded as numpy ndarrays), recursing into the elements
    if they contain matobjects.
    '''
    elem_list = []            
    for sub_elem in ndarray:
        if isinstance(sub_elem, sio.matlab.mio5_params.mat_struct):
            elem_list.append(_todict(sub_elem))
        elif isinstance(sub_elem,np.ndarray):
            elem_list.append(_tolist(sub_elem))
        else:
            elem_list.append(sub_elem)
    return elem_list
#
#def _todict(matobj):
#    '''
#    A recursive function which constructs from matobjects nested dictionaries
#    '''
#    dict = {}
#    for strg in matobj._fieldnames:
#        elem = matobj.__dict__[strg]
#        if isinstance(elem, np.ndarray):
#            i=1
#            for el in elem:
#                if isinstance(el, sio.matlab.mio5_params.mat_struct):
#                    dict[strg+'_'+str(i)] = _todict(el) 
#                    i+=1
#                else:
#                    dict[strg] = elem
#        elif isinstance(elem, sio.matlab.mio5_params.mat_struct):
#            dict[strg] = _todict(elem)
#        else:
#            dict[strg] = elem
#    return dict
#
# 

def stkshow(data):
    from pyqtgraph.Qt import QtCore, QtGui
    import pyqtgraph as pg
    import sys

    
    # determine if you need to start a Qt app. 
    # If running from Spyder, the answer is a no.
    # From cmd, yes. From nb, not sure actually.
    if not QtGui.QApplication.instance():
        app = QtGui.QApplication([])
    else:
        app = QtGui.QApplication.instance()
        
    ## Create window with ImageView widget
    win = QtGui.QMainWindow()
    win.resize(800,800)
    imv = pg.ImageView()
    win.setCentralWidget(imv)
    win.show()
    win.setWindowTitle('Fetching image stack...')
    
    
    
    
    resizeflg = 0;
    maxxysize = 800;
    maxdataxySc = np.floor(max(data.shape[0],data.shape[1])/maxxysize).astype('int')
    if maxdataxySc>1:
        resizeflg = 1;

       
    from numba import jit
    
    @jit(nopython = True)
    def DownScale(imgin): #use 2x downscaling for scrol speed   
       # imgout = downscale_local_mean(imgin,(Sc, Sc)
       imgout = (imgin[0::2,0::2]+imgin[1::2,0::2]+imgin[0::2,1::2]+imgin[1::2,1::2])/4
       return imgout
   
    
#    %different scales... might need anti aliasing - slow. stick to fast stuff for now
#    import math              
#    from numba import jit
#    
#    @jit
#    def mean_array(inp):
#        I, J = inp.shape
#    
#        mysum = 0
#        for i in range(I):
#            for j in range(J):
#                mysum += inp[i, j]
#        mysum = mysum/((I+1)*(J+1))
#        return mysum
#    
#    @jit
#    def AAandDownScale(imgin,Sc):    
#       # imgout = downscale_local_mean(imgin,(Sc, Sc))
#       imgout = np.zeros((np.ceil(imgin.shape/np.array([Sc,Sc]))).astype('int'),dtype = 'uint16')       
#       for ix in range(math.ceil(imgin.shape[0]/Sc)):
#           for iy in range(math.ceil(imgin.shape[1]/Sc)):
#               imgout[ix,iy] = mean_array(imgin[Sc*ix:Sc*(ix+1),Sc*(iy):Sc*(iy+1)])
#       return imgout
   
    
    if len(data.shape)==4:#RGB assume xytc
        if data.shape[3]==3 or data.shape[3]==4:
            if resizeflg:
                dataRs = np.zeros((np.ceil(data.shape/np.array([maxdataxySc,maxdataxySc,1,1]))).astype('int'),dtype = 'uint16')    
                for i in range(0,data.shape[2]):
                    for j in range(0,data.shape[3]):
                        dataRs[:,:,i,j] = DownScale(data[:,:,i,j])
                dataRs = dataRs.transpose((2,0,1,3))
            else:
                dataRs = data;
                dataRs = dataRs.transpose((2,0,1,3))
        else:
            sys.exit('color channel needs to be RGB or RGBA')
    elif len(data.shape)==3:
        if resizeflg:
            dataRs = np.zeros((np.ceil(data.shape/np.array([maxdataxySc,maxdataxySc,1]))).astype('int'),dtype = 'uint16')    
            for i in range(0,data.shape[2]):
                dataRs[:,:,i] = DownScale(data[:,:,i])
            dataRs = dataRs.transpose([2,0,1])
        else:
            dataRs = data;
            dataRs = dataRs.transpose([2,0,1])
                
    elif len(data.shape)==2:
        if resizeflg:
            dataRs = np.zeros((np.ceil(data.shape/np.array([maxdataxySc,maxdataxySc]))).astype('int'),dtype = 'uint16')
            dataRs = DownScale(data)
        else:
            dataRs = data;
    else:
        print('Data must be 2D image or 3D image stack')
    

    
    # Interpret image data as row-major instead of col-major
    pg.setConfigOptions(imageAxisOrder='row-major')
    

    win.setWindowTitle('Stack')
    
    ## Display the data and assign each frame a 
    imv.setImage(dataRs)#, xvals=np.linspace(1., dataRs.shape[0], dataRs.shape[0]))

    ##must return the window to keep it open
    return win
    ## Start Qt event loop unless running in interactive mode.
    if __name__ == '__main__':
        import sys
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtGui.QApplication.instance().exec_()
