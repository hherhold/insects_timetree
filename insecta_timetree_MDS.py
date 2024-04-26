#!/usr/bin/env python
from ete3 import Tree
from collections import defaultdict

#newick_file = "Insecta_order.nwk"
newick_file = "Insecta_family.nwk"
#newick_file = "Insecta_genus.nwk"

# make a data prefix string to use for the output files. This is the newick filename without the extension.
data_prefix = newick_file.split(".")[0]

t = Tree(newick_file,format = 1)

# bottom to top : use max
for n in t.traverse('postorder'):
    if n.is_leaf():
        n.height = 0
    else:
        n.height = max([c.height + c.dist for c in n.children])
# top to bottom to correct non-ultrametricity
for n in t.traverse('preorder'):
    if n.is_root():
        continue
    n.height = n.up.height - n.dist

def build_dist_dict( n , DD ):
    if n.is_leaf():
        DD[n.name][n.name] = n.height
        return [n.name],DD
    
    Cs = []
    for c in n.children:
        cs,DD = build_dist_dict( c , DD )
        Cs.append(cs)

    
    for i in range(len(Cs)-1):
        for j in range(i+1, len(Cs)):
            for c1 in Cs[i]:
                for c2 in Cs[j]:
                    DD[c1][c2] = n.height
                    DD[c2][c1] = n.height
                    
    C_total = []
    for cs in Cs:
        C_total += cs
    
    return C_total , DD
        
        
from collections import defaultdict

DD = defaultdict(dict)
_,DD = build_dist_dict( t , DD )

# Make a dataframe from the dictionary. The rows and columns are sorted. The diagonal is zero.
# Each cell is the distance between the row and column taxon.
import pandas as pd
D = pd.DataFrame(DD).sort_index().T.sort_index()

# Save the dataframe to a CSV file using the data prefix string.
D.to_csv(data_prefix + "_distances.csv")

import numpy as np
import scipy.spatial.distance as dist

## note to self : I checked this functions results versus R's cmdscale

# thanks to Marco Galardini and John Lees for this function
# source:https://pyseer.readthedocs.io/en/master/_modules/pyseer/cmdscale.html
def cmdscale(D , k , quiet=False):
    """Classical multidimensional scaling (MDS)

    Args:
        D (numpy.array)
            Symmetric distance matrix (n, n)

    Returns:
        Y (numpy.array)
            Configuration matrix (n, p). Each column represents a dimension. Only the
            p dimensions corresponding to positive eigenvalues of B are returned.
            Note that each dimension is only determined up to an overall sign,
            corresponding to a reflection.
        e (numpy.array)
            Eigenvalues of B (n, 1)
    """
    # Number of points
    n = len(D)

    # Centering matrix
    H = np.eye(n) - np.ones((n, n))/n

    # YY^T
    B = -H.dot(D**2).dot(H)/2

    # Diagonalize
    evals, evecs = np.linalg.eigh(B)

    # Sort by eigenvalue in descending order
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]

    # Compute the coordinates using positive-eigenvalued components only
    w, = np.where(evals > 0)
    L = np.diag(np.sqrt(evals[w]))
    V = evecs[:, w]
    Y = V.dot(L)

    if Y.shape[1]>=k:
        Y = Y[:, :k ]

    elif not quiet:
        print("Warning: MDS for {} dimensions, but only {} positive eigenvalues.".format(k,Y.shape[1]))


    return Y, evals[w]

D0 = D.copy()

# Convert the distance matrix to a numpy array.
D0 = D.values

D0[ np.diag_indices( D0.shape[0] ) ] = 0

Yfull,evals  = cmdscale(D , D0.shape[0] ) 

from sklearn.manifold import MDS


mds = MDS( n_components=2 , 
          eps=10**-6 , 
          max_iter = 100000 , dissimilarity = "precomputed" )
mds.fit( D0 )
mds.n_iter_

filename = data_prefix + ".mMDS.xy.csv"

pd.DataFrame( mds.embedding_  , columns = ['x','y'] , index = D.index).to_csv(filename)

mds3 = MDS( n_components=3 , 
          eps=10**-6 , 
          max_iter = 100000 , dissimilarity = "precomputed" )
mds3.fit( D0 )
mds3.n_iter_

mds3 = MDS( n_components=3 , 
          eps=10**-6 , 
          max_iter = 100000 , dissimilarity = "precomputed" )
mds3.fit( D0 )
mds3.n_iter_

# Save the MDS coordinates to a CSV file using the data prefix string.
filename = data_prefix + "_mds3.xyz.csv"
pd.DataFrame( mds3.embedding_  , columns = ['x','y','z'] , index = D.index).to_csv(filename)

