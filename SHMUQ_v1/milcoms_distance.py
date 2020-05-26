#
# compute distance between two milcoms files
# Jordan Fox 3/20

import numpy as np
import sys

def load_milcoms(filename):
    with open(filename) as fh:
        lines = np.array([float(x.strip()) for x in fh.readlines()])
        n = int(lines[0])
        eigvals = lines[1:(n+1)]
        eigvecs = lines[(n+1):]
        eigvecs = eigvecs.reshape((n,n))
        eigvecs = eigvecs.T
    return n, eigvals, eigvecs

def dist_metric(S,n):
    #d = np.linalg.norm(R - np.eye(subspace_dim))/np.linalg.norm(R)
    Sdiag,_ = np.linalg.eig(S)
    d = np.linalg.norm(Sdiag - np.ones(n))/n
    return d

subspace_dim = 5 # compute overlap matrix of first subspace_dim vectors

A_fn = sys.argv[1]
B_fn = sys.argv[2]

_, _, A = load_milcoms(A_fn)
_, _, B = load_milcoms(B_fn)

A = A[:,0:subspace_dim]
B = B[:,0:subspace_dim]

R = A.T @ B
print("A.T @ B =")
for r in R:
    print(r)
dist = dist_metric(R,subspace_dim)
print("distance = %3.3f" % dist)

make_plot = True
if make_plot:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.imshow(R)
    plt.colorbar()
    plt.savefig('overlap_matrix.pdf')

