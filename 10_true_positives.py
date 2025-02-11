import pickle as pkl
import networkx as nx
import random as rd
import pandas as pd



#Load ground truth data
file = open('./real_data/ground_truth_digraph.pkl','rb')
true_G = pkl.load(file)
file.close()

t_edges = [(u,v) for u,v in true_G.edges()]
total = len(t_edges)

#Load confusion matrix
data = pd.read_csv("./results/16s_graph_comparison.csv")


print("True graph has ",total," edges .")

#print("\t Best true positives:")
#print(data.sort_values('tp',ascending=False))

print("\t Best accuracy:")
print(data.sort_values('accuracy',ascending=False))

print("\t Best recall:")
print(data.sort_values('recall',ascending=False))

print("\t Best precision:")
print(data.sort_values('precision',ascending=False))

