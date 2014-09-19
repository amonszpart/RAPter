"""@package Primitive
This module provides an abstraction of the relationGraph using networkX

"""

import networkx as nx
import packages.primitive as primitive

class RelationGraph(object):
    def __init__(self,primArray, assignArray):
        self.G=nx.Graph()
        
        # First create the nodes
        for p in primArray:
            self.G.add_node(p.uid, w=0)
        
        # Then their relations
        for idx1, p1 in enumerate(primArray):
            for idx2, p2 in enumerate(primArray):
                if (idx2 > idx1):           
                    self.G.add_edge(p1.uid, p2.uid)
        
        # And finaly the node weights (number of samples)
        # assignArray[][0] = point id
        # assignArray[][0] = primitive uid
        for a in assignArray:
            self.G.node[a[1]]['w'] += 1
                    
        print "Number of primitives:  ",self.G.number_of_nodes()
        print "Number of connections: ",self.G.number_of_edges()
        
    def draw(self):
        nx.draw_spectral(self.G)
        #nx.draw_networkx_labels(self.G,pos=nx.spectral_layout(self.G))
