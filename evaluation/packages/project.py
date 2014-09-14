"""@package Project
This module defines an interface to read and use InputGen projects

See C++ InputGen project for more details on Projects
"""

import xml.etree.ElementTree as ET
import displacementKernels as kernels

class PyProject:
    """Main class of the module
    """

    def __init__(self, path):
        """Default constructor, taking as input a prj file."""
        print 'Loading project ' + path.name

        tree = ET.parse(path)
        root = tree.getroot()

        self.kernels = []

        #print kernels.generateDisplacementKernel(0)
        
        for groupNode in root:
            # extract displacement kernels
            if groupNode.tag == 'displacements':
                for kernelNode in groupNode:
                    if kernelNode.tag == 'kernel':
                        self.kernels.append( kernels.generateDisplacementKernel (kernelNode.attrib))

        print "Loaded kernels: "
        for k in self.kernels:
            print "   ", k
                

