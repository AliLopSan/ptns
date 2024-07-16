# Using real life data
import pickle as pkl
import pandas as pd
#our modules:
import tools as ts
import matrix as mx
from treebased import TB_Network


#Import species tree data
data = pd.read_csv("./real_data/interphylum_species_50.csv")
data['name_on_ncbi_tree'] = data['name_on_ncbi_tree'].str.replace("'","")
print (data.head())

#Create a dated species tree
S = ts.get_ultrametric_tree("./real_data/inter_phylum.phy",data)
print(S.to_newick())

file = open("./real_data/species_tree_ultra.pkl","wb")
pkl.dump(S,file)
file.close()

#Create character-state matrix
K0_list = ['K13628','K00561','K19310','K19115']
m = mx.generate_from_KEGG_v1(list(data['kegg_id']),K0_list)
file1 = open('./data/interphylum_matrix.pkl', 'wb')
pkl.dump(m,file1)
file1.close()

#Create a TB Network
N = TB_Network(S.root)
N.init_base_from_tralda(S,len(K0_list))

for leaf in N.leaves():
    leaf.chars = m[leaf.label]

ts.print_cs_matrix(N)

#Generate a basic Fitch Labeling
N.fitch_labeling()
fas = N.get_fas_by_state_change()

print("\t CHARACTER \t #fas \t fas")
for i in range(0,len(K0_list)):
    print("\t ",K0_list[i],"\t",len(fas[i]),"\t",fas[i])
