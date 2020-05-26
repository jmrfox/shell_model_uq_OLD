
# ShMUQ : Approximate Hessian module
# functions for computing approx. Hessian, covariance, MILCOMs, etc
# Fox 4/2020

import numpy as np

nl = "\n"

# FUNCTIONS

def truncate_data(data,method):
    if method is None:
        return data
    elif method == 'evenA':
        return data[(data[:,0]+data[:,1]) % 2 == 0,:]
    elif method == 'oddA':
        return data[(data[:,0]+data[:,1]) % 2 == 1,:]
    elif method == 'lowA':
        Amax = 32
        return data[(data[:,0]+data[:,1])<=Amax,:]
    elif method == 'ZNequal':
        return data[(data[:,0] == data[:,1]),:]
    elif method == 'gs':
        data_gs = []
        data_zno = data[:,[0,1,8]]
        data_zno_unique = np.unique(data_zno,axis=0)
        #print(data_zno_unique)
        #data_zno = [list(x) for x in data_zno]
        #data_zno_unique = [list(x) for x in data_zno_unique]
        for zno in data_zno_unique:
            #pdb.set_trace()
            subset = data[np.all(data_zno==zno,axis=1),:]
            #print(subset)
            data_gs.append(subset[np.argmin(subset[:,6]),:])
        return np.asarray(data_gs)
    elif method == 'lowJ':
        twoJmax = 8
        print("max J = %3.1f" % (max(data[:,2])/2))
        print("J max = %3.1f" % (twoJmax/2))
        return data[(data[:,2] <= twoJmax),:]
    elif method == 'highJ':
        twoJmin = 8
        print("max J = %3.1f" % (max(data[:,2])/2))
        print("J min = %3.1f" % (twoJmin/2))
        return data[(data[:,2] >= twoJmin),:]
    elif method == 'lowRes':
        Rmax = 1.
        return data[np.abs((data[:,6] - data[:,5])/de_notrunc) <= Rmax,:]
    elif method == 'highRes':
        Rmin = 1.
        return data[np.abs((data[:,6] - data[:,5])/de_notrunc) > Rmin,:]


def compute_jacobian(data,n_ops,theory_err):  # jacobian of observables wrt parameters
    exp_err_weighting = True
    n_data = data.shape[0]
    n_states = n_data//n_ops
    states = data[:n_states,[0,1,2,3,4]] #z n j t n
    de = np.sqrt(data[:n_states,7]**2 + (theory_err*np.ones(n_states))**2)
    znjtno = data[:,[0,1,2,3,4,8]] #z n j t n op#
    jac = np.zeros([n_ops,n_states])
    vec_i = np.zeros(n_states)
    for i in range(n_ops):
        for j,st in enumerate(states):
            st_oi = np.append(st,i+1)
            current_state = data[(znjtno==st_oi).all(axis=1),:][0]
            if exp_err_weighting:
                de = np.sqrt(current_state[7]**2 + theory_err**2)
                vec_i[j] = current_state[9]/de
            elif not exp_err_weighting:
                vec_i[j] = current_state[9]
        jac[i,:] = vec_i
    return jac

def compute_hessian(J):  # hessian approx of chi-squared wrt parameters
    hess = J @ J.T
    return hess

def write_matrix(M,out_fn,dim):
    with open(out_fn,'w') as out_fh:
        out_fh.write(str(dim)+nl)
        for i in range(dim):
            for j in range(dim):
                out_fh.write("%i\t%i\t%f\n" % (i+1,j+1,M[i,j]))
        print('File written: '+out_fn)

def compute_obs_cov(J,C):
    Cobs = J.T @ C @ J
    return Cobs

def load_matrix(fn_in,dim,header=False):
    #loads a matrix in format: i j Aij
    # returns a numpy array
    A = np.zeros([dim,dim])
    with open(fn_in,'r') as fh_in:
        counter = 0
        for line in fh_in:
            counter += 1
            if header and (counter == 1):
                continue
            else:
                ls = line.split()
                A[int(ls[0])-1,int(ls[1])-1] = float(ls[2])
    return A

def compute_milcoms(A,fn_out):
    w, v = np.linalg.eig(A)   # v = eigenvecs in columns
    idx = w.argsort()[::-1]
    w = w[idx]
    v = v[:,idx]
    n = len(w)
    with open(fn_out,'w') as fh_out:
        fh_out.write(str(n) + nl)
        for wi in w:
            #print(wi)
            fh_out.write(str(wi) + nl)
        for i in range(n): #column idx
            for j in range(n): #row idx
                fh_out.write(str(v[j,i]) + nl)
    print("File written: %s" % fn_out)
