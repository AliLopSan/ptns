import random as rd
import asymmetree.treeevolve as te

from treebased import TB_Node,TB_Network

#Choose an arbitrary number of species and characters
n_species = rd.randint(3,25)
n_chars   = rd.randint(1,25)
print("\t Generate a random PTN on ",n_species," species and ",n_chars," characters.")

#Generate a random species tree with n_species
S = te.species_tree_n(n_species,planted=False)
print("\t Species Tree: ",S.to_newick())

#Initialize TB_Network object
N = TB_Network(S.root)
N.init_base_from_tralda(S,n_chars)

#Generate character-state matrix for leaves
character_state = dict()
for species in S.leaves():
    temp = []
    for i in range(0,n_chars):
        coin_flip = rd.randint(0,1)
        if coin_flip == 0:
            temp.append(False)
        else:
            temp.append(True)
    character_state[species.label] = list(temp)

for leaf in N.leaves():
    leaf.chars = character_state[leaf.label]


#Generate a Fitch labeling for the inner nodes
N.fitch_labeling()

#Run Greedy algorithm
fa_list = N.get_fas_by_state_change()
N.greedy_completion(fa_list,S)
print("\t Greedy completion added: ",len(N.transfers))
N.print_transfer_hw_info()
