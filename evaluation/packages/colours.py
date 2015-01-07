"""@package Colours
This module provides the colourmaps used in globOpt to display primitives 
according to their gid

"""

import packages.primitive as primitive
import packages.orderedSet as orderedSet

class Colours(object):

    def __init__(self):
        self.colListMedium = ['#F15A60', '#7AC36A', '#5A9BD4', '#FAA75B', '#9E67AB', '#CE7058', '#D77FB4', '#F1ADCB', '#B2A377']
        self.colListDark = ['#F15A60', '#7AC367', '#5A9B15', '#FAA75B', '#9E67AB', '#CE7058', '#D77FB4', '#F1ADCB', '#B2A377']
    
    """
    Compute the colourmap associating one colour per group id. Also output the
    masks associating the node idx for each group id (can be used directly as
    filter in networkX display funtions)
    """
    def getDIDColourMap(self, primArray):
        ids = orderedSet.OrderedSet()
        gfilter = {}
        for p in primArray:
            ids.add(p.did)
            
            if p.did not in gfilter:
                gfilter[p.did] = []
            gfilter[p.did].append(p.uid)
        
        cmap    = {}
        nbCol   = len(self.colListMedium)
        
        for idx, did in enumerate(ids):
            print idx, idx%nbCol, self.colListMedium[idx%nbCol]
            cmap[did] = self.colListMedium[idx%nbCol]
        
        return cmap, gfilter
