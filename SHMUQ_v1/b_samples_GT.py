#
#   print out N samples (.int) given variances of MILCOMs
#   Fox 2019
#

import numpy as np
from sa_mod import *
import glob, os
from collect_strengths import *
import pickle
from trans_mod import get_trans_data

class nucleus:
    def __init__(self,element,Zv,Nv,twoJz):
        self.element = element
        self.Zv = Zv
        self.Nv = Nv
        self.A = core_mass + self.Zv + self.Nv
        self.name = self.element + str(self.A)
        self.scaling = ['1',str(core_mass+2),str(self.A),'0.3']
        self.twoJz = twoJz

### SET PARAMETERS HERE

n_samples = 1000
use_exact_interaction = False
do_bigstick_runs = True
parallel = True
collect_only = False

interaction_name = 'usdb'
core_mass = 16
sps = "sd"
parent = nucleus('Si',6,10,0)
daughter = nucleus('P',7,9,0)
Tshift = 2
T2int = "T2sd"
#opme_name = "GTsd_0.78"
op_type = "GT"
opme_name = "GTsd"
isospin_mirror = False
if parent.Zv==daughter.Nv and daughter.Zv==parent.Nv:
    isospin_mirror = True

diag_opt = "ld"
nkeep_n = 20
nkeep_d = 20

####
nthreads=8
nranks = 2

gtstrength_cmd = "./gtstrength.x"
if parallel:
    include_frag = True
    os.environ['OMP_NUM_THREADS']=str(nthreads)
    #bigstick_cmd = 'mpirun -n '+str(nranks)+' bigstick-mpi-omp.x'
    bigstick_cmd = './bigstick-openmp.x'
elif not parallel:
    include_frag = False
    bigstick_cmd = './bigstick.x'

fn_out_pkl = parent.name+"_"+opme_name+".pkl"

"""
def get_trans_data(case_name):
    # transition data. trans[<case>] is a list of dictionaries with keys 2Ji, ni, 2Jf, nf , B
    # 2Ji and 2Jf are initial and final 2J
    # ni and nf are 1,2,3,... corresponding to the first, second, third, ... state with that J-value

    trans = {}
    trans[case_name] = [{'2Ji':0,'ni':1,'2Jf':2,'nf':1,'B':0.0},\
            ]

    return trans
"""

def make_sample_interactions(n_samples):

    print('Writing %i sample interaction files...' % n_samples)
    hess_approx_evals = np.loadtxt('milcoms_approx.dat',\
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

def make_bigstick_inputs(int_name_list,nuc,opt):

    print('Writing bigstick inputs...')
    #int_files = []
    #for int_fn in glob.glob("usdb_rand*.int"):
    #    int_files.append(int_fn[:-4])

    #opt = "d"
    ZvNv = " ".join([str(nuc.Zv),str(nuc.Nv)])
    twoJz = str(nuc.twoJz)
    frag = "0"
    #diag_opt = "ld"
    if opt=='d':
        nkeep = nkeep_d
    elif opt=='n':
        nkeep = nkeep_n
    lanit = str(20*nkeep)
    nkeep = str(nkeep)
    lanopt = " ".join([nkeep,lanit])
    #lanopt = nkeep

    run_name_list = []
    for int_name in int_name_list:
        if opt=='n':
            run_name = nuc.name+int_name
        elif opt=='d':
            if Tshift>0:
                run_name = nuc.name+'_Tsh_'+int_name
            else:
                run_name = nuc.name+int_name
        scaling = nuc.scaling
        in_list = [opt,run_name,sps,ZvNv,twoJz]
        if include_frag:
            in_list.append(frag)
        in_list = in_list + [int_name," ".join(scaling)]
        if Tshift>0 and opt=='d':
            tshift_scaling = [-1*Tshift,-1*Tshift,0,0]
            in_list = in_list + [T2int,' '.join([str(x) for x in tshift_scaling])]
        in_list = in_list + ["end",diag_opt,lanopt]

        outfn = "in."+run_name+".b"
        with open(outfn,'w') as outfh:
             outfh.write("\n".join(in_list)+"\n")
        run_name_list.append(run_name)
    return run_name_list

def run_bigstick(run_name_list):
    for run_name in run_name_list:
        print('Running bigstick for '+run_name)
        cmd = " ".join([bigstick_cmd,"<",'in.'+run_name+'.b'])
        subprocess.call(cmd,shell=True)

def make_gtstrength_inputs(int_name_list):
    print('Writing gtstrength inputs...')

    #res_files = []
    #for fn in glob.glob(parent.name+"*.res"):
    #    res_files.append(fn[:-4])

    #if len(res_files) == 0:
    #    print("No matching files")

    #hbarc = 197.3
    #hbaromega = 41.0 * A ** (-1./3.)
    #mass = 0.5 * (938.272 + 939.565)
    #scaling = hbarc / np.sqrt(hbaromega * mass)
    #scaling = A ** (1/3)  # b^2
    ## close approx :  b = A^(1/6)
    ## for Cr48, b = 1.0


    run_name_list = []
    for int_name in int_name_list:
        if isospin_mirror:
            par_name = parent.name+int_name
            dau_name = par_name
            den_name = par_name
        else:
            par_name = parent.name+int_name
            dau_name = daughter.name+int_name
            if Tshift>0:
                den_name = daughter.name + '_Tsh_' + int_name
            else:
                den_name = dau_name

        par_ZN = ' '.join([str(x) for x in [parent.Zv,parent.Nv]])
        dau_ZN = ' '.join([str(x) for x in [daughter.Zv,daughter.Nv]])

        run_name = par_name+opme_name
        outfn = "in." + run_name +".gstr"
        in_list = [opme_name,"1.0",par_name,"0",dau_name,"0",\
                den_name,str(Tshift),"n",run_name,par_ZN,dau_ZN]
        with open(outfn,'w') as outfh:
            outfh.write("\n".join(in_list)+"\n")
        run_name_list.append(run_name)
    return run_name_list

def run_gtstrength(run_name_list):
    for run_name in run_name_list:
        input_fn = 'in.'+run_name+'.gstr'
        print('Running gtstrength for '+input_fn)
        cmd = " ".join([gtstrength_cmd,"<",input_fn])
        subprocess.call(cmd,shell=True)

def parse_str(int_name,case_name,trans_data):
    fn = case_name+int_name+opme_name+'.str'  #filename convention
    print("Parsing "+fn+"...")
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
                if line_num < 3:
                    continue
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
        results_dict[int_name] = get_trans_data(parent.name,op_type)
        results_dict = parse_str(int_name,parent.name,results_dict)
    print(results_dict)
    f = open(fn_out_pkl,"wb")
    pickle.dump(results_dict,f)
    f.close()

if __name__=="__main__":
    if collect_only:
        int_name_list=[]
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
        if Tshift>0:
                run_name_list_par = make_bigstick_inputs(int_name_list,parent,'n')
                run_name_list_dau = make_bigstick_inputs(int_name_list,daughter,'n')
                run_name_list_den = make_bigstick_inputs(int_name_list,daughter,'d')
                run_bigstick(run_name_list_par)
                run_bigstick(run_name_list_dau)
                run_bigstick(run_name_list_den)
        elif isospin_mirror:
            run_name_list_par = make_bigstick_inputs(int_name_list,parent,'d')
            #run_name_list_dau = make_bigstick_inputs(int_name_list,daughter,'n')
            #run_name_list_den = make_bigstick_inputs(int_name_list,daughter,'d')
            run_bigstick(run_name_list_par)
            #run_bigstick(run_name_list_dau)
            #run_bigstick(run_name_list_den)
        else:
            run_name_list_par = make_bigstick_inputs(int_name_list,parent,'n')
            #run_name_list_dau = make_bigstick_inputs(int_name_list,daughter,'n')
            run_name_list_den = make_bigstick_inputs(int_name_list,daughter,'d')
            run_bigstick(run_name_list_par)
            #run_bigstick(run_name_list_dau)
            run_bigstick(run_name_list_den)

    run_name_list_g = make_gtstrength_inputs(int_name_list)
    run_gtstrength(run_name_list_g)

    collect_strengths(int_name_list)
    print("Done")




