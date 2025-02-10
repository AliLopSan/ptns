#..........................................................
#  Testing Sankoff's labeling with different weights
#..........................................................
import pickle as pkl
import pandas as pd
import sys
import os
from treebased import TB_Node,TB_Network


# Creates a dataframe for KOs
def create_ko_df(name,ko_list,fa_dict):
    fa_set = []
    fa_count = []
    for i in range(0,len(ko_list)):
        fa_list = [node.label for node in fa_dict[i]]
        fa_list = sorted(fa_list)
        fa_set.append(fa_list)
        fa_count.append(len(fa_list))

    df = pd.DataFrame()
    df['characters'] = ko_list
    df['fa_length'] = fa_count
    df['fa_list'] = fa_set
    print(df.head())
    df.to_csv(name,index=False)


#Import Species Tree
f = open("./real_data/species_tree_ultra.pkl","rb")
S = pkl.load(f)
f.close()

#Import KO list
KO_df = pd.read_csv("./real_data/ARG_related_KOs.csv")
print(KO_df.head())
print(KO_df.shape)
KOs_list = list(KO_df['RF_KO'])
matrix = [KOs_list[i:i+45] for i in range(0,len(KOs_list),45)]


#Sankoff labeling weights
penalization = [(1.0,1.0),(0.5,1.0),(1.0,0.5),(0.25,1),(1,0.25)]
penaliz_type = ["equal","loss_half","hgt_half","loss_quarter","hgt_quarter"]

for j in range(0,len(penalization)):
    folder_name = "./real_data/sankoff/penalizations/" + penaliz_type[j] + "/"
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    loss_cost,fa_cost = penalization[j]
    
    for i in range(0,len(matrix)):
        chunk = matrix[i]
        #import character-state matrix
        mname = "./real_data/KO_matrices/interphylum_matrix_chunk_" + str(i) + ".pkl" 
        file = open(mname,"rb")
        m = pkl.load(file)
        file.close()
    
        #Initialize tree-based network
        N = TB_Network(S.root)
        N.init_base_from_tralda(S,len(chunk))

        for leaf in N.leaves():
            leaf.chars = m[leaf.label]
            
        #Generate a sankoff labeling for the inner nodes
        N.sankoff_labeling(loss_cost,fa_cost)

        #Calculate first appearance nodes
        fas = N.get_fas_by_state_change()
        dname = folder_name + "fa_info_dataset_"+str(i)+"_SAN.csv" 
        create_ko_df(dname,chunk,fas)


#Merge penalizations
for j in range(0,len(penalization)):
    folder_name = "./real_data/sankoff/penalizations/" + penaliz_type[j] + "/"
    sankoff_data = pd.DataFrame()
    for i in range(0,4):
        name = folder_name + 'fa_info_dataset_' + str(i) + '_SAN.csv'
        indata = pd.read_csv(name)
        sankoff_data = pd.concat([sankoff_data, indata])
    fname = folder_name + "FA_sankoff_"+ penaliz_type[j] +".csv"
    sankoff_data.to_csv(fname,index=False)

