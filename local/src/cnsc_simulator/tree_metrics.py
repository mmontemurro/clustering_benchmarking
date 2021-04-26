#!/usr/bin/env python3

import sys, os, argparse
import numpy as np
from gen_tree import gen_tree 

######################################################
#                                                    #
# Utility functions developed by Marilisa Montemurro #
# to compute metrics on a tree                       #
#                                                    #
######################################################
"""
class MyNode(Node):
    def __init__(self, name, parent=None):
        Node.__init__(self, name, parent)
        self.id = 0
        self.name = name
        self.parent=parent
        self.tuple=[]
        self.is_dead=False
        self.edge_length = 0
        # alelle length for each chromosome, root has the same as reference
        self.cn=[]
        self.chrlen=[]
        self.ref=[]
        self.snvs=[]
        self.corres=[]
        self.cn_summary={}
        self.cn_detail=[]
        self.parentID = -1
        self.subTree=-1
        self.subTree_l2=-1
    def getTuple(self):
        return self.tuple
    def setDead(self):
        self.is_dead=True
    def getID(self):
        return self.id

"""
  
# Function to check if there is a path from root 
# to the given node. It also populates 
# 'arr' with the given path 
def getPath(tree, rarr, n):
      
    # push the node's id in 'arr' 
    rarr.append(n.id)
  
    # if it is the root node 
    # return true 
    if n.parent.id == -1:
        return True
      
    # else recur on the parent node
    if getPath(tree, rarr, tree[n.parent.id]):
        return True
      
    return False


def distanceBetweenNodes(tree, n1, n2):
    # vector to store the path of 
    # first node n1 to root 
    path1 = []
  
    # vector to store the path of 
    # second node n2 to root 
    path2 = []
    getPath(tree, path1, tree[n1])
    getPath(tree, path2, tree[n2])

    # Get intersection point
    i, j = 0, 0
    intersection=-1
    path1.reverse() #from root to n1
    path2.reverse() #from root to n2

    while(i != len(path1) or j != len(path2)):
  
        # Keep moving forward until no intersection 
        # is found
        if (i == j and path1[i] == path2[j]):
            i += 1
            j += 1
        else:
            intersection = j - 1
            break
  
    # Save path
    path = []
    for i in range(len(path1) - 1, intersection - 1, -1):
        path.append(path1[i])
    for j in range(intersection + 1, len(path2)):
        path.append(path2[j])

    return len(path) - 1

# Driver program 
if __name__=='__main__':
  
    parser = argparse.ArgumentParser("Utility functions to compute metrics on trees")
    parser.add_argument("-i", "--input", required=True, help="Path of tree npy structure")
    parser.add_argument("-d", "--distance", action="store_true", help="Compute distance between two nodes")
    parser.add_argument("-n", "--nodes_id", required= '--distance' in sys.argv, help="Identifiers of the nodes", nargs=2, type=int)

    args = parser.parse_args()

    tree = np.load(args.input, allow_pickle=True)
    d = distanceBetweenNodes(tree, args.nodes_id[0], args.nodes_id[1])
    print("Path between {} and {} lenght = {}".format(args.nodes_id[0], args.nodes_id[1], d))

    