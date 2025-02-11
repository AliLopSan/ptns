import pickle as pkl
import networkx as nx
import random as rd
import pandas as pd


def build_infered_graph(hw_data,t_threshold,sdata,S):
    G = nx.DiGraph()
    kegg_2_name = dict(zip(sdata['kegg_id'].to_list(),sdata['species'].to_list()))
    name_2_index  = dict()
    name_2_leaves = dict()
    leaves = []
    i = 0
    
    #Creates nodes
    for node in S.leaves():
        name_2_index[node.label] = i
        leaves.append(node.label)
        i = i + 1
        G.add_node(node.label)
        G.nodes[node.label]['species'] = node.label
        
    #Creates a dict that maps nodes to leaves in underlying subtree
    for node in S.postorder():
        if node.is_leaf():
            name_2_leaves[node.label] = [node.label]
        else:
            leaves = []
            for child in node.children:
                leaves = leaves + name_2_leaves[child.label]
            name_2_leaves[node.label] = leaves
    data = pd.read_csv(hw_data)
    x = list(data['a'])
    y = list(data['b'])
    w = list(data['weight'])
    for i in range(0,len(w)):
        if w[i] >= t_threshold:
            x_leaves = name_2_leaves[x[i]]
            y_leaves = name_2_leaves[y[i]]

            for x_leaf in x_leaves:
                for y_leaf in y_leaves:
                    G.add_edge(x_leaf,y_leaf)
                    G.add_edge(y_leaf,x_leaf)
    
    return(G)
    
def contingency_table(true_graph,graph):
    t_nodes = [n for n in true_graph.nodes()]
    t_edges = [(u,v) for u,v in true_graph.edges()]
    i_nodes = [graph.nodes[n]['species'] for n in graph.nodes()]
    i_edges = [(graph.nodes[u]['species'],graph.nodes[v]['species']) for u,v in graph.edges()]

    if (len(t_nodes) != len(i_nodes) or set(t_nodes) != set(i_nodes)):
        raise ValueError("compared graphs must have the same vertex sets")

    tp,fp,fn = 0,0,0

    for u,v in i_edges:
        if (u,v) in t_edges:
            tp +=1
        else:
            fp +=1
    for u,v in t_edges:
        if (u,v) not in i_edges:
            fn += 1
    tn = (len(i_nodes) * (len(i_nodes)-1) - (tp + fp + fn))
    
    return tp,tn,fp,fn


def performance(true_graph,graph):
    tp,tn,fp,fn = contingency_table(true_graph,graph)
    accuracy = (tp + tn) / (tp + tn + fp + fn) if tp + tn + fp + fn > 0 else float('nan')
    precision = tp / (tp + fp) if tp + fp > 0 else float('nan')
    recall = tp / (tp + fn) if tp + fn > 0 else float('nan')
    f1 = (precision*recall)/(precision + recall) if precision + recall > 0 else float('nan')
    
    return (graph.order(), graph.size(),
            tp, tn, fp, fn,
            accuracy, precision, recall,f1)




#Load ground truth data
file = open('./real_data/ground_truth_digraph.pkl','rb')
true_G = pkl.load(file)

#Loading our data
method = ['sankoff','genesis']
penaliz_type = ["equal","hgt_half","hgt_quarter"]
t_thresholds = [5,9,18,27,36,45]

f1 = open("./real_data/species_tree_ultra.pkl","rb")
S  = pkl.load(f1)
f1.close()

sdata = pd.read_csv("./real_data/interphylum_species_50.csv")

m_id = []
p_id = []
th_list = []
tp_list = []
tn_list = []
fp_list = []
fn_list = []
acc_list = []
precision_list = []
recall_list = []
f1_list = []


for m in method:
    for p in penaliz_type:
        data_name = "./results/"+m+"_"+ p +"_hw_info.csv"
        for t in t_thresholds:
            G = build_infered_graph(data_name,t,sdata,S)
            order,size,tp,tn,fp,fn,accuracy,precision,recall,f1 = performance(true_G,G)
            print("\t For dataset: ",m," ",p," ")
            print("\ttp: ",tp," tn: ",tn," fp: ",fp,"fn:",fn)
            print("\taccuracy: ",accuracy,"\t precision",precision,"\t recall",recall)
            m_id.append(m)
            p_id.append(p)
            th_list.append(t)
            tp_list.append(tp)
            tn_list.append(tn)
            fp_list.append(fp)
            fn_list.append(fn)
            acc_list.append(accuracy)
            precision_list.append(precision)
            recall_list.append(recall)
            f1_list.append(f1)
                
#Add the Fitch dataset               
data_name = "./results/fitch_hw_info.csv"

for t in t_thresholds:
    G = build_infered_graph(data_name,t,sdata,S)
    order,size,tp,tn,fp,fn,accuracy,precision,recall,f1 = performance(true_G,G)
    m_id.append("fitch")
    p_id.append("nope")
    th_list.append(t)
    tp_list.append(tp)
    tn_list.append(tn)
    fp_list.append(fp)
    fn_list.append(fn)
    acc_list.append(accuracy)
    precision_list.append(precision)
    recall_list.append(recall)
    f1_list.append(f1)

results = pd.DataFrame()
results['labeling method'] = m_id
results['penalization type'] = p_id
results['transfer_threshold'] = th_list
results['tp'] = tp_list
results['tn'] = tn_list
results['fp'] = fp_list
results['fn'] = fn_list
results['accuracy'] = acc_list
results['precision'] = precision_list
results['recall'] = recall_list
results['F1'] = f1_list



results.to_csv("./results/16s_graph_comparison.csv",index=False)

print(results)
