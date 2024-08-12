import matplotlib.pyplot as plt
import seaborn as sns
import pickle as pkl
import numpy as np
import pandas as pd
from treebased import TB_Node,TB_Network

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

def create_ko_df(ko_list,fa_dict):
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
    df.to_csv('./real_data/fa_info_dataset.csv',index=False)

    

def plot_KO_hist(ko_list,fa_dict):
    leaves_count = []
    inner_count  = []

    for i in range(0,len(ko_list)):
        leaves = 0
        inner  = 0
        for node in fa_dict[i]:
            if node.is_leaf():
                leaves+=1
            else:
                inner+=1
        leaves_count.append(leaves)
        inner_count.append(inner)

    sns.set()
    #For legend
    fa_type = ['leaves','internal nodes']
    x_pos = np.arange(len(ko_list))
    bar_width = 0.35
    
    #Two bars should be created
    plt.bar(x_pos, leaves_count, color=(0.5,0.1,0.5,0.6))
    plt.bar(x_pos, inner_count, color='blue',bottom=leaves_count)
    
    # Add title and axis names
    plt.title(r'First appearances per character', fontsize=20)
    plt.xlabel(r'characters',fontsize=15)
    plt.ylabel(r'number of first appearance nodes')
    
    
    plt.xticks(x_pos, ko_list,rotation='vertical',fontsize=10)
    plt.legend(fa_type)
    plt.axhline(y = 2, color = 'g', linestyle = '--',lw=2.5) 
 
    # Show graph
    plt.savefig('FA_distribution_counts.pdf',format='pdf',bbox_inches='tight')
    plt.show()
    plt.close()

#Import Species Tree
f = open("./real_data/species_tree_ultra.pkl","rb")
S = pkl.load(f)
f.close()

#Import Character-State Matrix
file = open("./real_data/interphylum_matrix_23_KOs.pkl","rb")
m = pkl.load(file)
file.close()


#The list of KOs related to antibiotic resistance
f1_K0s = ['K18220','K02257','K01610','K04068','K00627','K00241','K01679','K13628','K01669','K03980','K06886','K07305','K01589','K00561','K07483','K17836','K02227','K18214','K19310','K19115','K02274','K03737','K00850']


#Initialize TB_Network object
N = TB_Network(S.root)
N.init_base_from_tralda(S,len(f1_K0s))


for leaf in N.leaves():
    leaf.chars = m[leaf.label]


#Generate a Fitch labeling for the inner nodes
N.fitch_labeling()

#Calculate first appearance nodes
fas = N.get_fas_by_state_change()

#Plot & Calculate first-appearance information
for i in range(0,len(f1_K0s)):
    print("\t ",f1_K0s[i],"\t",len(fas[i]),"fas according to Fitch")

plot_KO_hist(f1_K0s,fas)
create_ko_df(f1_K0s,fas)
