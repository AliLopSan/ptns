import pickle as pkl
import numpy as np
import pandas as pd
import os
from treebased import TB_Node,TB_Network

def genesis_labeling(TB,loss_cost,fa_cost):
    def _base_case_value(dp_table,v,label,gain,char):
        if v.chars[char] == False:
            if label == 0 and gain == 0:
                return(0)
            else:
                return(1000000000) #not inf but something very big
        else:
            if label == 1:
                return(0)
            else:
                return(1000000000) #same

    def _inner_case_value(dp_table,v,current_label,current_gain,M):
        #Special case: V[v,0,1]
        # This case fixes one child as origin and leaves the rest to fate
        if current_label == 0 and current_gain == 1:
            child_sum = []
            for child in v.children:
                origin_penalty = []
                for l in [0,1]:
                    origin_penalty.append(dp_table[child][l,1] + M[current_label,l])
                brother_sum = 0
                for brother in v.children:
                    if brother != child:
                        label_cost = []
                        for l in [0,1]:
                            label_cost.append(dp_table[brother][l,0] + M[current_label,l])
                        brother_sum = brother_sum + min(label_cost)
                child_sum.append(min(origin_penalty) + brother_sum)
            return(min(child_sum))
        else:
            cost = 0
            for child in v.children:
                label_costs = []
                for l in [1,0]:
                    label_costs.append(dp_table[child][l,0] + M[current_label,l])
                cost = cost + min(label_costs)
            return cost

    #Backtracking
    def _find_best_child_origin(dp_table,v,M):
        argmin = []
        children = list(v.children)
        for child in v.children:
            origin_penalty = []
            for l in [0,1]:
                origin_penalty.append(dp_table[child][l,1] + M[0,l])
            brother_sum = 0
            for brother in v.children:
                if brother != child:
                    label_cost = []
                    for l in [0,1]:
                        label_cost.append(dp_table[brother][l,0] + M[0,l])
                    brother_sum = brother_sum + min(label_cost)
            argmin.append(min(origin_penalty) + brother_sum)
        
        index_min = np.argmin(argmin)
        return(children[index_min])


    #cost matrix
    M = np.array([[0.0,fa_cost],[loss_cost,0.0]])
    dp_table = dict()

    # Main algo
    for char in range(0,len(TB.root.chars)):
        dp_table = dict()

        #main algo
        for v in TB.postorder():

            dp_table[v] = dict()
            for label in [0,1]:
                for gain in [0,1]:
                    if v.is_leaf():
                        dp_table[v][label,gain] = _base_case_value(dp_table,v,label,gain,char)
                    else:
                        dp_table[v][label,gain] = _inner_case_value(dp_table,v,label,gain,M)

        #backtracking
        for v in TB.preorder():
            if v == N.root:
                if dp_table[v][1,1] <= dp_table[v][0,1]:
                    v.chars[char] = True
                    origin = N.root
                else:
                    v.chars[char] = False
                    #choose the origin
                    origin = _find_best_child_origin(dp_table,v,M)
           
            else:
                if v == origin:
                    if dp_table[v][1,1] <= dp_table[v][0,1]:
                        v.chars[char] = True
                    else:
                        v.chars[char] = False
                        #choose the origin
                        if not v.is_leaf():
                            origin = _find_best_child_origin(dp_table,v,M)
                else:
                    if dp_table[v][1,0] <= dp_table[v][0,0]:
                        v.chars[char] = True
                    else:
                        v.chars[char] = False



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
    folder_name = "./real_data/genesis/penalizations/" + penaliz_type[j] + "/"
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
        genesis_labeling(N,loss_cost,fa_cost)

        #Calculate first appearance nodes
        fas = N.get_fas_by_state_change()
        dname = folder_name + "fa_info_dataset_"+str(i)+"_GEN.csv" 
        create_ko_df(dname,chunk,fas)


#Merge penalizations
for j in range(0,len(penalization)):
    folder_name = "./real_data/genesis/penalizations/" + penaliz_type[j] + "/"
    sankoff_data = pd.DataFrame()
    for i in range(0,4):
        name = folder_name + 'fa_info_dataset_' + str(i) + '_GEN.csv'
        indata = pd.read_csv(name)
        sankoff_data = pd.concat([sankoff_data, indata])
    fname = folder_name + "FA_genesis_"+ penaliz_type[j] +".csv"
    sankoff_data.to_csv(fname,index=False)
