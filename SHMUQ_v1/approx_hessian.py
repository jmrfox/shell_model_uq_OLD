#!/usr/bin/python


# ShMUQ: Hessian calculation
# computes approximate hessian, param. covariance, MILCOMs

import numpy as np
#import multiprocessing as mp
import os
#import itertools
#import pdb
import argparse
from approx_hessian_mod import *

if __name__=="__main__":

    ##theory_err = 0.158545 # set reduced X^2 = 1 for USDB

    nl = "\n"

    parser = argparse.ArgumentParser(description="Computes approximate Hessian matrix, parameter covariance matrix, and PCA-parameters from input expectation values")
    parser.add_argument('exp_values',type=str,nargs=1,help="file of expectation values")
    parser.add_argument('output_tag',type=str,nargs=1,help="Tag (string) to denote outputs")
    parser.add_argument('-tu','--theory_unc',type=float,nargs=1,help="Theory error in MeV",default=0.1)
    parser.add_argument('-n','--n_params',type=int,nargs=1,help="Number of parameters in interaction",default=66)
    parser.add_argument('-tr','--trunc',nargs=1,help="Data truncation method. See hessian_approx_mod for details.",default=None)

    args = parser.parse_args()

    want_obs_cov = False
    want_milcoms = True

    fn_in = args.exp_values[0]
    data = np.loadtxt(fn_in,delimiter="\t",skiprows=1)
    ###  data = Z(0)   N(1)   J(2)   T(3)   n(4)   Esm(5)    Eexp(6)     dEexp(7)   op#(8)     <o>(9)
    trunc_method = args.trunc
    if trunc_method is not None:
        data = truncate_data(data,trunc_method[0])

    output_tag = args.output_tag[0]
    fn_out_hes = 'hessian_approx_'+output_tag+'.dat'
    fn_out_cov = 'covariance_approx_'+output_tag+'.dat'
    fn_out_mil = 'milcoms_approx_'+output_tag+'.dat'

    theory_unc = args.theory_unc
    print("ATTENTION: the fixed theoretical uncertainty is %f" % theory_unc)
    n_ops = args.n_params

    J = compute_jacobian(data,n_ops,theory_unc)
    A = compute_hessian(J)
    write_matrix(A,fn_out_hes,n_ops)
    C = np.linalg.inv(A)
    write_matrix(C,fn_out_cov,n_ops)

    if want_milcoms:
        compute_milcoms(A,fn_out_mil)

    if want_obs_cov:
        Cobs = compute_obs_cov(J,C)
        write_matrix(Cobs,'obs_cov_approx_'+output_tag+'.dat')



