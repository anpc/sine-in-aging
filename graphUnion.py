import networkx as nx
import pickle
import matplotlib.pyplot as plt

filename1 = 'graphPart1'  # graph
filename2 = 'graphPart2'  # graph
filename3 = 'graphPart3'  # graph
filename4 = 'graphPart4'  # graph
filename5 = 'graphPart5'  # graph
filename6 = 'graphPart6'  # graph
filename7 = 'graphPart7'  # graph


infile1 = open(filename1, 'rb')
G1 = pickle.load(infile1)

infile2 = open(filename2, 'rb')
G2 = pickle.load(infile2)

infile3 = open(filename3, 'rb')
G3 = pickle.load(infile3)

infile4 = open(filename4, 'rb')
G4 = pickle.load(infile4)

infile5 = open(filename5, 'rb')
G5 = pickle.load(infile5)

infile6 = open(filename6, 'rb')
G6 = pickle.load(infile6)

infile7 = open(filename7, 'rb')
G7 = pickle.load(infile7)

listGraph= [G1, G2,G3, G4,G5,G6,G7]

graphUnion = nx.compose_all(listGraph)

cliqueCounter=0
notCliqueCounter=0
maxclique=0
maxNoclique=0
hist = {}
histnotClique = {}

comps=nx.connected_components(graphUnion)

for comp in comps:
        singel_comp = graphUnion.subgraph(comp).copy()
        nodeNum = singel_comp.number_of_nodes()
        edgeNum = singel_comp.number_of_edges()
        if (edgeNum ==  ( nodeNum * ( nodeNum - 1 )) / 2 ) :
#                hist[nodeNum] = hist.get(nodeNum, 0) + 1
                cliqueCounter += 1
                if (nodeNum > maxclique):
                    maxclique= nodeNum
        else :
                histnotClique[nodeNum] = histnotClique.get(nodeNum, 0) + 1
                notCliqueCounter +=1
                if (nodeNum == 23):
                    hist[nx.graph_clique_number(singel_comp)] = hist.get(nx.graph_clique_number(singel_comp), 0) + 1
                if (nodeNum > maxNoclique):
                    maxNoclique= nodeNum
                    biggestnoC= singel_comp

#pickle.dump(listGraph, outfile)
#outfile.close()
print (hist)
print(histnotClique)
print ('biggest clique', maxclique)
print ('biggest not clique', maxNoclique)
#print ('biggest clique in naxnoclique', nx.node_clique_number(biggestnoC))
print('size of qlique: ', cliqueCounter)
print('size of not qlique: ',notCliqueCounter)
infile1.close()
infile2.close()
infile3.close()
infile4.close()
infile5.close()
infile6.close()
infile7.close()

