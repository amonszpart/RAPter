"""@package DisplacementKernels
This module defines an interface to read and use InputGen projects

See C++ InputGen project for more details on Projects
"""
from numpy import vectorize
import math
import numpy as np
from scipy import stats


class DisplacementKernel(stats.rv_continuous):
    """Python representation of the C++ class AbstractDisplacementKernel
    """
    def __init__(self, typeId = None, sname = None, enabled = None, *args, **kwargs):   
        super(DisplacementKernel, self).__init__(*args, **kwargs)
        
    def getName(self):
        return self._ctor_param['sname']

    def getType(self):
        return self._ctor_param['typeId']

    def isEnabled(self):
        return self._ctor_param['enabled']
        
    def _stats(self):
        return 0., 0., 0., 0.

    def __str__(self):
        return "%s displacement kernel (enabled=%s)" % (self._ctor_param['sname'], self._ctor_param['enabled'])
        

class UniformRandomDisplacementKernel(DisplacementKernel):
    """Python representation of the C++ class UniformRandomDisplacementKernel
    """
    def __init__(self, paramList, enabled, *args, **kwargs):
        super(UniformRandomDisplacementKernel, self).__init__(*args, **kwargs)        
        self._ctor_param['paramList'] = paramList
        self._ctor_param['sname']   = "Random (Uniform)"
        self._ctor_param['typeId']  = 0
        self._ctor_param['enabled'] = enabled
        
    def _pdf(self, x):
        rangeMin = self._ctor_param['paramList'][0]
        rangeMax = self._ctor_param['paramList'][1]
        pdfvalue = 1. / (rangeMax - rangeMin)
        
        return np.where((x<rangeMin) or (x>rangeMax), 0., )
    

class NormalRandomDisplacementKernel(DisplacementKernel):
    """Python representation of the C++ class NormalRandomDisplacementKernel
    """
    def __init__(self, paramList, enabled, *args, **kwargs):
        super(NormalRandomDisplacementKernel, self).__init__(*args, **kwargs)         
     
        self._ctor_param['paramList'] = paramList
        self._ctor_param['sname']   = "Random (Normal)"
        self._ctor_param['typeId']  = 1
        self._ctor_param['enabled'] = enabled
        #self._ctor_param['momtype'] = 0
    
    def _pdf(self, x):
        mean  = self.mean()
        stdev = self.stdev()
        
        #print 'did you call me ?'

        #http://www.cplusplus.com/reference/random/normal_distribution/
        return 1./(stdev* math.sqrt(2.*math.pi)) * np.exp(-pow((x-mean),2)/(2.*stdev*stdev))
        
    def mean(self): 
        return self._ctor_param['paramList'][1]
        
    def stdev(self): 
        return self._ctor_param['paramList'][0]
            

class BiasDisplacementKernel(DisplacementKernel):
    """Python representation of the C++ class BiasDisplacementKernel
    """
    def __init__(self, bias, enabled = True, *args, **kwargs):
        super(BiasDisplacementKernel, self).__init__("Bias", 2, enabled, *args, **kwargs)         

        self.bias  = bias

    # The PDF of this function is a Dirac pic. 
    def _pdf(self, xarray):
        print "Warning: Bias PDF is an infinite Dirac that will not been displayed"
        return 0

"""Factory to generate displacement kernels from its typeId
    \param paramArray Array generated when parsing the xml file
"""
def generateDisplacementKernel(paramArray):
    def UniformRandomDisplacementParam(paramArray):
        return float(paramArray['distributionMin']), float(paramArray['distributionMax'])
    def NormalRandomDisplacementParam(paramArray):
        return float(paramArray['distributionStdDev']), float(paramArray['distributionMean'])
    def BiasDisplacementParam(paramArray):
        return float(paramArray['bias'])

    def isKernelEnabled(paramArray):
        return paramArray['enabled'] != '0'

    factory = { '0': UniformRandomDisplacementKernel,
        '1': NormalRandomDisplacementKernel,
        '2': BiasDisplacementKernel }
    paramfactory = { '0': UniformRandomDisplacementParam,
        '1': NormalRandomDisplacementParam ,
        '2': BiasDisplacementParam }

    typeId = paramArray['typeId']

    return factory[typeId](paramfactory[typeId](paramArray), isKernelEnabled(paramArray))

