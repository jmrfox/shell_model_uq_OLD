#
#   print out N samples (.int) given variances of MILCOMs
#   Fox 2019
#

import numpy as np
from sa_mod import *
import glob, os
#from collect_strengths import *
import pickle
import os
from hocalc_mod import oscb
from trans_mod import get_trans_data

n_samples = 1000
#milcoms_file ='milcoms_augmented.dat'
milcoms_file ='milcoms_approx.dat'

interaction_name = 'usdb'  # leave off .int (does not need to be sorted)

Zv = 5
Nv = 5
nuc_name = "Al"
core_mass = 16
sps = "sd"
A=core_mass+Zv+Nv
twoJz = str((Zv+Nv)%2)
#twoJz = str(2)

diag_opt = 'ex'
nkeep = "10"
lanit = "400"

op_type = "E2"
if op_type == "M1":
    opme_name = "M1sd"
    scaling = 1.0
elif op_type == "E2":
    #opme_name = "E2sd_1.57_0.45"    # Mg26 optimized
    #opme_name = "E2sd_1.29_0.49"
    #opme_name = "E2sd_1.35_0.65"
    opme_name = "E2sd_1.36_0.45"
    #opme_name = "E2sd_1.5_0.5"
    #opme_name = "E2sd_1.57_0.45"
    #opme_name = "E2sd_1.6_0.6"
    #opme_name = "E2sd_1.65_0.35"
    #scaling = A ** (1/3)  # b^2
    scaling = oscb(A)**2
    #scaling = 1.802**2      #Mg26 from Brown-Chung-Wildenthal 1980

genstrength_cmd = "./genstrength.x"
parallel = True
nranks = 1
nthreads = 8
if parallel:
    diag_opt = 'ld'
    include_frag = True
    os.environ['OMP_NUM_THREADS']=str(nthreads)
    #bigstick_cmd = 'mpirun -n '+str(nranks)+' bigstick-mpi-omp.x'
    bigstick_cmd = './bigstick-openmp.x'
elif not parallel:
    include_frag = False
    bigstick_cmd = './bigstick.x'

use_exact_interaction = False
collect_only = False
do_bigstick_runs = True    #for debug only

case_name = nuc_name+str(A)
fn_out_pkl = case_name+"_"+opme_name+".pkl"

#def get_trans_data():
#    # transition data. trans[<case>] is a list of dictionaries with keys 2Ji, ni, 2Jf, nf , B
#    # 2Ji and 2Jf are initial and final 2J
#    # ni and nf are 1,2,3,... corresponding to the first, second, third, ... state with that J-value
#
#    trans = {}
#    trans[case_name] = [{'2Ji': 4  ,'ni': 1   ,'2Jf': 0    ,'nf': 1 ,'B':0.0},
#                        {'2Ji': 4  ,'ni': 2   ,'2Jf': 0    ,'nf': 1 ,'B':0.0},
#                        {'2Ji': 0  ,'ni': 2   ,'2Jf': 4    ,'nf': 1 ,'B':0.0},
#                        {'2Ji': 4  ,'ni': 4   ,'2Jf': 0    ,'nf': 1 ,'B':0.0},
#                        ]
#
#    return trans

def make_sample_interactions(n_samples):
    print('Writing %i sample interaction files...' % n_samples)
    hess_approx_evals = np.loadtxt(milcoms_file,\
            skiprows=1,delimiter="\t")[:66]
    param_variance = 1/hess_approx_evals

#usdbmil = np.loadtxt('usdbmil.vec', skiprows=1)
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

def make_bigstick_inputs(int_name_list):

    print('Writing bigstick inputs...')
    #int_files = []
    #for int_fn in glob.glob("usdb_rand*.int"):
    #    int_files.append(int_fn[:-4])

    opt = "d"
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

def run_bigstick(input_fn):
    print('Running bigstick for '+input_fn)
    cmd = " ".join([bigstick_cmd,"<",input_fn])
    subprocess.call(cmd,shell=True)

def make_genstrength_inputs():
    print('Writing genstrength inputs...')

    res_files = []
    for fn in glob.glob(case_name+"*.res"):
        res_files.append(fn[:-4])

    if len(res_files) == 0:
        print("No matching files")

    #hbarc = 197.3
    #hbaromega = 41.0 * A ** (-1./3.)
    #mass = 0.5 * (938.272 + 939.565)
    #scaling = hbarc / np.sqrt(hbaromega * mass)
    #scaling = A ** (1/3)  # b^2
    ## close approx :  b = A^(1/6)
    ## for Cr48, b = 1.0
    out_fn_list = []
    for fn in res_files:
        print(fn)

        outfn = "in."+fn+opme_name+".gstr"
        in_list = [opme_name,str(scaling),fn,"0",fn,"0",fn,"n",fn+opme_name]
        with open(outfn,'w') as outfh:
            outfh.write("\n".join(in_list)+"\n")
        out_fn_list.append(outfn)
    return out_fn_list

def run_genstrength(input_fn):
    print('Running genstrength for '+input_fn)
    cmd = " ".join([genstrength_cmd,"<",input_fn])
    errcode = subprocess.call(cmd,shell=True)
    print('genstrength error code = %i' % errcode)
    return errcode

def parse_str(int_name,case_name,trans_data):
    fn = case_name+int_name+opme_name+'.str'  #filename convention
    trans_list = trans_data[int_name][case_name]
    for io,o in enumerate(trans_list):
        line_num = 0
        with open(fn,'r') as fh:
            print('transition',o)
            counteri = o['ni'] - 1
            counterf = o['nf'] - 1
            found_parent = False
            for line in fh:
                line_num += 1

                ls = line.split()
                if ('parent' in ls) and (float(ls[2])*2 == o['2Ji']):
                    if (counteri==0):
                        found_parent = True
                        print("found parent :",line)
                    elif (counteri>0):
                        counteri = counteri - 1
                if found_parent and ('daughter' in ls) and (float(ls[2])*2 == o['2Jf']):
                    if (counterf==0):
                        trans_data[int_name][case_name][io]['B'] = float(ls[4])
                        print("found daughter :",line)
                        break
                    elif (counterf>0):
                        counterf = counterf - 1
    return trans_data

def collect_strengths(int_name_list):
    results_dict = {}    # results_dict[<interaction>] = trans_data updated with B-values
    for int_name in int_name_list:
        results_dict[int_name] = get_trans_data(case_name,op_type)
        results_dict = parse_str(int_name,case_name,results_dict)
    print(results_dict)
    f = open(fn_out_pkl,"wb")
    pickle.dump(results_dict,f)
    f.close()

if __name__=="__main__":
    if collect_only:
        int_name_list = []
        for i in range(n_samples):
            int_name = interaction_name+"_rand"+str(i).zfill(4)
            int_name_list.append(int_name)
        collect_strengths(int_name_list)
        exit('Done!')

    if use_exact_interaction:
        int_name_list = [interaction_name]
    else:
        int_name_list = make_sample_interactions(n_samples)

    if do_bigstick_runs:
        input_fn_list_b = make_bigstick_inputs(int_name_list)
        for input_fn in input_fn_list_b:
            run_bigstick(input_fn)

    input_fn_list_g = make_genstrength_inputs()
    for input_fn in input_fn_list_g:
        _ = run_genstrength(input_fn)

    collect_strengths(int_name_list)




