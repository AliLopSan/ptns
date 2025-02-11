import pandas as pd
import pickle as pkl
import numpy as np
import seaborn as sns
import sys
import matplotlib.pyplot as plt
from mne.viz.circle import _plot_connectivity_circle

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

def plot_fa_network_v2(sdata,data,S,name,t_threshold):
    def _create_phyl_tuples(sdata,name_2_index,leaves):
        colors = dict()
        #Dictionary for phylum color
        phyl_color = dict()
        phyl_color['Firmicutes'] = '#E52B50'
        phyl_color['Bacteroidetes'] = '#318CE7'
        phyl_color['Actinobacteria'] = '#8DB600'
        phyl_color['Proteobacteria'] = '#8A2BE2'

        KEGG = list(sdata['kegg_id'])
        PHY  = list(sdata['phylum'])

        for i in range(0,len(KEGG)):
            if PHY[i] in phyl_color.keys():
                colors[name_2_index[KEGG[i]]] = phyl_color[PHY[i]]
            else:
                colors[name_2_index[KEGG[i]]] = '#D3D3D3'
        temp = [colors[n] for n in range(0,len(leaves))]
                
        return temp

    kegg_2_name = dict(zip(sdata['kegg_id'].to_list(),sdata['species'].to_list()))
    name_2_index  = dict()
    name_2_leaves = dict()
    leaves = []
    i = 0
    for node in S.leaves():
        name_2_index[node.label] = i
        leaves.append(node.label)
        i = i + 1

    full_names = [ kegg_2_name[node] for node in leaves]

    for node in S.postorder():
        if node.is_leaf():
            name_2_leaves[node.label] = [node.label]
        else:
            leaves = []
            for child in node.children:
                leaves = leaves + name_2_leaves[child.label]
            name_2_leaves[node.label] = leaves

    list_of_lists = []

    for i in range(0,len(leaves)):
        row = []
        for i in range(0,len(leaves)):
            row.append(0)
        list_of_lists.append(row)
    
    x = list(data['a'])
    y = list(data['b'])
    w = list(data['weight'])

    for i in range(0,len(w)):
        if w[i] >= t_threshold:
            x_leaves = name_2_leaves[x[i]]
            y_leaves = name_2_leaves[y[i]]

            for x_leaf in x_leaves:
                x_pos = name_2_index[x_leaf]
                for y_leaf in y_leaves:
                    y_pos = name_2_index[y_leaf]
                    list_of_lists[x_pos][y_pos] = w[i]
    
    results = np.array(list_of_lists)
    con = np.where(results > 1, results, np.nan)
    colors  = _create_phyl_tuples(sdata,name_2_index,leaves)
    fig,axes = _plot_connectivity_circle(con, full_names,node_colors=colors, node_edgecolor=colors,
                                         colormap='viridis',facecolor='white',textcolor='black',fontsize_colorbar=20)
    fig.savefig(name,format='pdf',bbox_inches='tight')
    plt.show()
    plt.close()


f1 = open("./real_data/species_tree_ultra.pkl","rb")
S  = pkl.load(f1)
f1.close()


sdata = pd.read_csv("./real_data/interphylum_species_50.csv")


#Possible transfer highway thresholds
t_thresholds = [5,9,18,27,36,40]


#Start with Fitch dataset               
data_name = "./results/fitch_hw_info.csv"
data = pd.read_csv(data_name)
for t in t_thresholds:
        res_name = "./results/plots/pairwise_connectiivity/fitch_pairs_"+"t_"+str(t)+".pdf"
        plot_fa_network_v2(sdata,data,S,res_name,t)

#Parameters that make sense for the rest of the methods
method = ['sankoff','genesis']
penaliz_type = ["equal","hgt_half","hgt_quarter"]
sdata = pd.read_csv("./real_data/interphylum_species_50.csv")


for p in penaliz_type:
    data_name = "./results/sankoff_"+ p +"_hw_info.csv"
    data = pd.read_csv(data_name)
    for t in t_thresholds:
        res_name = "./results/plots/pairwise_connectiivity/sankoff_"+p+"_pairs_"+"t_"+str(t)+".pdf"
        plot_fa_network_v2(sdata,data,S,res_name,t)

for p in penaliz_type:
    data_name = "./results/genesis_"+ p +"_hw_info.csv"
    data = pd.read_csv(data_name)
    for t in t_thresholds:
        res_name = "./results/plots/pairwise_connectiivity/genesis_"+p+"_pairs_"+"t_"+str(t)+".pdf"
        plot_fa_network_v2(sdata,data,S,res_name,t)
