from __future__ import print_function

import networkx as nx
import packages.relationGraph as relgraph
import packages.primitive as primitive


def process(primArray1, assignArray1, 
            primArray2, assignArray2,
            primitiveCorres,
            primitiveCorresId,
            angles  = [],
            graph1  = None,
            graph2  = None,
            verbose = False):
            
    def mprint(str, *args):
        if verbose: print(str % args)

    ################################################################################
    ## Build relation graphs
    if graph1 == None:
        graph1 = relgraph.RelationGraph(primArray1, assignArray1)
    if graph2 == None:
        graph2 = relgraph.RelationGraph(primArray2, assignArray2)
    
    useAngles = len(angles) != 0

    for e in graph2.G.edges_iter(data=True):
        e[-1]['matched']=0
    for e in graph1.G.edges_iter(data=True):
        e[-1]['matched']=0

    for p in primArray1:
        p_node = graph1.G.node[p.uid]
        if not primitiveCorres.has_key(p):
            mprint ("Gt Primitive not matched (%i,%i)",p.uid,p.did)
        else:
            matched_p = primitiveCorres[p]
            matched_p_node = graph2.G.node[matched_p.uid]
            
            # iterate over all relations and check they have a counterpart in the estimated scene
            # cUid is the uid of the connected component
            # cUid can be used to access the connected component using
            # mprint (cUid, primitiveCorresId[cUid])
            for idx, cUid in enumerate(graph1.G.edge[p.uid]):
                # now we are looking for the connection starting from matched_p et going to primitiveCorresId[cUid]
                # matched_p.uid, primitiveCorresId[cUid]
                # 
                # if we find it, we increment the matched field of the edges, and move to the next one
                #mprint ((cUid, primitiveCorresId[cUid]))
                if cUid in primitiveCorresId:
                    for idx2, matched_cUid in enumerate(graph2.G.edge[matched_p.uid]):                
                        #mprint (("  ",matched_cUid, primitiveCorresId[cUid]))
                        if matched_cUid == primitiveCorresId[cUid]:
                            #mprint ("  match found !")
                            graph1.G.edge[p.uid][cUid]['matched'] += 1
                            graph2.G.edge[matched_p.uid][matched_cUid]['matched'] += 1
                            break
                else:
                    mprint("relation.py:54") 

            
    def checkEdges(graph, doPrint):
        correct=0
        error = False
        for e in graph.G.edges_iter(data=True):
            count = e[-1]['matched']
            if (count == 2):correct+=1
            elif count == 0: 
                if doPrint: mprint ("Missed edge detected...")
            else: error = True;
            
        return correct, error

    gtcorrect, gtError = checkEdges(graph1, True)
    correct_it1, error_it1 = checkEdges(graph2, False)

    if gtError or error_it1: 
        mprint ("Error occurred, invalid number of matches. ABORT")
        return -1., -1.
        
    if gtcorrect != correct_it1:
        mprint ("Error: non-symmetric detection")
        return -1., -1.
        
    mprint("Graph1, nbedges = %i", graph1.G.number_of_edges())
    mprint("Graph2, nbedges = %i", graph2.G.number_of_edges())
        
    precision = float(gtcorrect)/float(graph2.G.number_of_edges())
    recall    = float(gtcorrect)/float(graph1.G.number_of_edges())

    mprint ("precision = %f",precision)
    mprint ("recall    = %f",recall)
    
    return precision, recall
