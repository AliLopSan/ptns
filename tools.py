import pickle as pkl
import pandas as pd
from alive_progress import alive_bar
from asymmetree.tools.PhyloTreeTools import *
from tralda.datastructures.Tree import Tree
from treebased import TB_Node,TB_Network

# Functions to include:
# The count types of trianges should be an auxiliary function 
# measure fp,tn,tp, etc.
# visualization functions

#Print character-state matrix
def print_cs_matrix(N):
    print("\t Character-state Matrix")
    print("\t taxa \t\t characters")
    for leaf in N.leaves():
        print("\t ",leaf, "\t",leaf.chars)

#functions to be used with the matrix module
# parsing formats from other data types
#.........................................................
#    GET SPECIES TREE
#..........................................................
# INPUT: PHYLLIP format newick string generated by ncbi taxonomy common tree viewer
#       data which is a pandas df that maps ncbi id to keggid
# OUTPUT:  tralda tree from .phy file with a random ultrametric timing

def get_ultrametric_tree(phy,data):
    #Convert newick to tralda tree
    f = open(phy,"r")
    newick = ''
    for r1 in f:
        line = r1.replace('\n','')
        newick = newick + line
    tree = parse_newick(newick)
    #convert tralda tree into networkx with keggid labels
    name_2_keggid = dict(zip(data['name_on_ncbi_tree'],data['kegg_id']))
    
    for v in tree.preorder():
        if v.is_leaf():
            #name = v.label[0:v.label.index(':')]
            name = v.label.replace("'","")
            v.label = name_2_keggid[name]

    random_ultrametric_timing(tree,True,True)
    return tree
