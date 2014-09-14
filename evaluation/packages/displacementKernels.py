"""@package DisplacementKernels
This module defines an interface to read and use InputGen projects

See C++ InputGen project for more details on Projects
"""
from decimal import Decimal


class DisplacementKernel(object):
    """Python representation of the C++ class AbstractDisplacementKernel
    """
    def __init__(self, name, typeId, enabled):
        self.name    = name
        self.typeId  = typeId
        self.enabled = enabled

    def getName(self):
        return self.name

    def getType(self):
        return self.typeId

    def isEnabled(self):
        return self.enabled

    def __str__(self):
        return "%s displacement kernel (enabled=%s)" % (self.name, self.enabled)

class UniformRandomDisplacementKernel(DisplacementKernel):
    """Python representation of the C++ class UniformRandomDisplacementKernel
    """
    def __init__(self, paramList, enabled):
        super(UniformRandomDisplacementKernel, self).__init__("Random (Uniform)", 0, enabled)

        self.rangeMin = paramList[0]
        self.rangeMax = paramList[1]

class NormalRandomDisplacementKernel(DisplacementKernel):
    """Python representation of the C++ class NormalRandomDisplacementKernel
    """
    def __init__(self, paramList, enabled):
        super(NormalRandomDisplacementKernel, self).__init__("Random (Normal)", 1, enabled)         

        self.mean  = paramList[0]
        self.stdev = paramList[1]

class BiasDisplacementKernel(DisplacementKernel):
    """Python representation of the C++ class BiasDisplacementKernel
    """
    def __init__(self, bias, enabled):
        super(BiasDisplacementKernel, self).__init__("Bias", 2, enabled)         

        self.bias  = bias

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

