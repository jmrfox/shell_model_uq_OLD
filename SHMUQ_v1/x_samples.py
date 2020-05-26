#
#   compute sample interactions, and expectation values
#   Fox 2019
#

import numpy as np
from sa_mod import *
import glob, os
#from collect_strengths import *
import pickle
import os

n_samples = 5000

Zv = 10
Nv = 10
nuc_name = "Ar"
core_mass = 16
sps = "sd"
A=core_mass+Zv+Nv
twoJz = str((Zv+Nv)%2)
#twoJz = str(2)

diag_opt = 'ex'
nkeep = "10"
lanit = "200"

opme_name = 'sdls1bis'  # operator to compute exp val
interaction_name = 'usdb'
fn_milcoms = 'approx_milcoms.dat'

parallel = True
nranks = 1
nthreads = 16
if parallel:
    diag_opt = 'ld'
    include_frag = True
    os.environ['OMP_NUM_THREADS']=str(nthreads)
    bigstick_cmd = 'mpirun -n '+str(nranks)+' bigstick-mpi-omp.x'
elif not parallel:
    include_frag = False
    bigstick_cmd = './bigstick.x'

use_exact_interaction = False
collect_only = False
do_bigstick_runs = True    #for debug only

case_name = nuc_name+str(A)
fn_out_pkl = case_name+"_"+opme_name+".pkl"


def make_sample_interactions(n_samples):

    print('Writing %i sample interaction files...' % n_samples)
    hess_approx_evals = np.loadtxt(fn_milcoms,\
            skiprows=1,delimiter="\t")[:66]
    param_variance = 1/hess_approx_evals

    samples = np.random.multivariate_normal(mean = np.zeros(66),\
            cov = np.diag(param_variance), size = n_samples)

    int_name_list = []
    for i,s in enumerate(samples):
        int_name = interaction_name+"_rand"+str(i).zfill(4)
        pert = np.stack((np.arange(1,67),s),axis=-1)
        perturb_milcom(int_name,pert)
        int_name_list.append(int_name)
    return int_name_list

#os.chdir("/mydir")

def make_spectrum_inputs(int_name_list):

    print('Writing bigstick inputs...')
    opt = "n"
    ZvNv = " ".join([str(Zv),str(Nv)])
    #twoJz = str((Zv+Nv)%2)
    frag = "0"
    #diag_opt = "ld"
    #diag_opt = "ex"
    #nkeep = "20"
    #lanit = "400"
    lanopt = " ".join([nkeep,lanit])
    #lanopt = nkeep

    out_fn_list = []
    for int_name in int_name_list:
        outfn = "in."+case_name+int_name+".b"
        A = core_mass + Zv + Nv
        scaling = ["1.0",str(core_mass+2),str(A),"0.333"]
        if include_frag:
            in_list = [opt,case_name+int_name,sps,ZvNv,twoJz,\
                    frag,int_name," ".join(scaling),"end",diag_opt,lanopt]
        else:
            in_list = [opt,case_name+int_name,sps,ZvNv,twoJz,\
                    int_name," ".join(scaling),"end",diag_opt,lanopt]

        with open(outfn,'w') as outfh:
             outfh.write("\n".join(in_list)+"\n")
        out_fn_list.append(outfn)
    return out_fn_list

def make_expect_inputs(int_name_list):

    print('Writing expectation value inputs...')

    opt = "x"

    out_fn_list = []
    for int_name in int_name_list:
        wfn = case_name+int_name
        outfn = "inx."+case_name+int_name+".b"
        A = core_mass + Zv + Nv
        scaling = ["1.0"]*4
        in_list = [opt,wfn,wfn+opme_name,\
                opme_name," ".join(scaling),"end"]

        with open(outfn,'w') as outfh:
             outfh.write("\n".join(in_list)+"\n")
        out_fn_list.append(outfn)
    return out_fn_list

def run_bigstick(input_fn):
    print('Running bigstick for '+input_fn)
    cmd = " ".join([bigstick_cmd,"<",input_fn])
    subprocess.call(cmd,shell=True)

def read_x(int_name,case_name):
    fn = case_name+int_name+opme_name+'.res'  #filename convention
    spectrum = getspectrum(fn)
    print(spectrum)
    return float(spectrum[0][-2])


def collect_x(int_name_list):
    results = []    # results list
    for int_name in int_name_list:
        results.append(read_x(int_name,case_name))
    print(results)
    f = open(fn_out_pkl,"wb")
    pickle.dump(results,f)
    f.close()

if __name__=="__main__":
    if collect_only:
        int_name_list = []
        for i in range(n_samples):
            int_name = interaction_name+"_rand"+str(i).zfill(4)
            int_name_list.append(int_name)
        collect_x(int_name_list)
        exit('Done!')

    if use_exact_interaction:
        int_name_list = [interaction_name]
    else:
        int_name_list = make_sample_interactions(n_samples)

    if do_bigstick_runs:
        input_fn_list_s = make_spectrum_inputs(int_name_list)
        for input_fn in input_fn_list_s:
            run_bigstick(input_fn)
        input_fn_list_x = make_expect_inputs(int_name_list)
        for input_fn in input_fn_list_x:
            run_bigstick(input_fn)

    collect_x(int_name_list)




