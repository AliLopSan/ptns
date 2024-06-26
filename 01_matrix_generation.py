import asymmetree.treeevolve as te
import matrix as mx
import random as rd

from treebased import TB_Network

def print_cs_matrix(N):
    print("\t Character-state Matrix")
    print("\t taxa \t\t characters")
    for leaf in N.leaves():
        print("\t ",leaf, "\t",leaf.chars)


#Choose an arbitrary number of species and characters.........................
n_species = rd.randint(3,10)
n_chars   = rd.randint(1,10)
print("\t Generate a PTN on ",n_species," species and ",n_chars," characters.")

#Generate a species tree with asymmetree.....................................
S = te.species_tree_n(n_species,planted=False)
print("\t Species Tree: ",S.to_newick())
leaves = [s.label for s in S.leaves()]

#Generate a PTN with random matrix............................................
n1 = TB_Network(S.root)
n1.init_base_from_tralda(S,n_chars)

matrix = mx.generate_random_matrix(S,n_chars)

for leaf in n1.leaves():
    leaf.chars = matrix[leaf.label]
print("\t Randomly ")
print_cs_matrix(n1)

#Generate a PTN with simple HGT model .......................................
#hgt_rate = rd.uniform(0,1)
hgt_rate = 0.5
loss_rate = 0.9
#loss_rate = rd.uniform(0,1)
params = {'HGT_type':"simple",
          "dup_rate":0.0,
          "hgt_rate":hgt_rate,
          "loss_rate":loss_rate}
name = "simple"
n2 = TB_Network(S.root)
n2.init_base_from_tralda(S,n_chars)
matrix,highway = mx.generate_from_asymm(S,n2,n_chars,params,name)

for leaf in n2.leaves():
    leaf.chars = matrix[leaf.label]

print("\t Using simple HGT model:")
print("\t hgt rate:",hgt_rate,"\t loss rate:",loss_rate)
print_cs_matrix(n2)
highway.print_transfer_hw_info()

#Generate a PTN with inverse model   ........................................
params["HGT_type"] = "inverse"
name = "inverse"
n3 = TB_Network(S.root)
n3.init_base_from_tralda(S,n_chars)
matrix,highway = mx.generate_from_asymm(S,n3,n_chars,params,name)

for leaf in n3.leaves():
    leaf.chars = matrix[leaf.label]

print("\t Using inverse bias HGT model:")
print("\t hgt rate:",hgt_rate,"\t loss rate:",loss_rate)
print_cs_matrix(n3)
highway.print_transfer_hw_info()
