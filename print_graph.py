import networkx as nx
import pickle
import matplotlib.pyplot as plt


filename = 'graph'  # graph

infile = open(filename, 'rb')
G2 = pickle.load(infile)

infile.close()

nx.draw(G2)
plt.show()
