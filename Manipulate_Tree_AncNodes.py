#modifies a tree nodes and leaves
from ete3 import Tree

### Not Run: to add species names in the tree. 
# species_dict = {}
# with open('/Users/timeisen/Dropbox (Personal)/KuriyanLab/ImmuneCellSignaling/Analysis/Analysis20240322/species/idmapping_2024_05_06.tsv', 'r') as f:
#    next(f) #skip header
#    for line in f:
#       organism = line.split("\t")[5]
#       organism = organism.replace("(", "<")
#       organism = organism.replace(")", ">")
#       species_dict[line.split("\t")[0]] = organism

tree = Tree('output_tree_TEST.tre')
# print(tree.write(format=8))
for node in tree.traverse('preorder'):
   if not node.is_leaf():
      node.name = str(int(node.support)) + "_" + str(int(node.support) + 243)
      # node.name = str(int(node.support) + 243)



for node in tree.traverse('preorder'):
   if node.is_leaf():
      # print(node.support)
      node.name = "_".join(node.name.split("_")[1:]).split("/")[0]
      # node.name = species_dict[node.name]

# Dictionary to store extracted subtrees
extracted_subtrees = {}

# Traverse the tree and identify the internal nodes with the desired IDs
for node in tree.traverse():
   if not node.is_leaf():
        node_id = node.name
        if node_id == '143_386': #TXK
            extracted_subtrees[node_id] = node.copy()

combined_tree = Tree()
for node_id, subtree_root in extracted_subtrees.items():
    combined_tree.add_child(subtree_root)


print(combined_tree.write(format = 8, format_root_node=False))
