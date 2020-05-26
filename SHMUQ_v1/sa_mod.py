# Sensitivity Analysis module for ShMUQ

import numpy as np
import os, sys
import subprocess
import time
import pdb
import glob
import operator

verb = True
#verb = False

# parameters
nkeep=40
nkeep_repeat=20
liter=800

fn_milcoms = "approx" #filename for milcoms file without .milcom ext
fn_milcoms_vector = "usdbmil.vec"
fn_int_vec_sorted = 'usdb_sorted.vec'
interaction_name = 'usdb'

#theory_err = 0.1 # from USD paper
theory_err = 0.158545 # set reduced X^2 = 1

# BINARIES
bigstick_parallel = "./bigstick-openmp.x"
bigstick_serial = "./bigstick.x"
anamil_cmd = "./anamil.x"
standard_cmd = "./standardorderint.x"

#brown_data_fn = 'BrownData.csv'
#exp_data_fn = 'Ne21Data.csv'

#
#   DATA
#

def get_data(data_fn,max_unc=2.0,abs_energies=False,swap_jt=True):
    # data columns = A(0) , Z(1) , N(2) , 2t(3) , 2j(4) , n(5) , energy(6) , error(7)
    exp_data = np.genfromtxt(data_fn, delimiter=',',skip_header=1)

    if max_unc is not None:   # delete data with dE >= max_unc in MeV
        exp_data = exp_data[exp_data[:,7]<max_unc]

    if not abs_energies:
        # if abs_energies is TRUE then all energies are absolute
        # else, assume that g.s. energy is absolute and the rest are excitations
        for i,line in enumerate(exp_data):  #absolute energies
            if line[6]<0:     # if find g.s.
                E0 = line[6]
            elif line[6]>0:    # add g.s. energy to excitations to get absolute
                exp_data[i,6] = line[6] + E0

    if swap_jt:   # if you need to swap j and t columns
        exp_data[:,[3,4]] = exp_data[:,[4,3]]   # swap 2t 2j  ->  2j 2t

    # get interaction  original values
    int_vec = np.genfromtxt(fn_int_vec_sorted,skip_header=1)
    # exp_data = A(0), Z(1) , N(2) , 2j(3) , 2t(4) , n(5) , energy(6) , error(7)
    # where energy(6) is absolute energy
    # THIS IS THE STANDARD COLUMN ORDER
    return exp_data,int_vec

#
#   LITTLE TOOLS
#

def times_so_far(ls):
    out = [0]*len(ls)
    for i in range(len(ls)):
        out[i] = ls[:i+1].count(ls[i])
    return out

def round_to_half(x):
    return round(2*x)/2

def getspectrum(filename):
    if verb: print('Getting spectrum from '+filename)
    outs = []
    with open(filename,"rt") as fptr_in:
        contents = fptr_in.read()
        lines = contents.split("\n")
        istate = 1
        for line in lines:
            try:
                if (len(line.split())==5 or len(line.split())==6) and int(line.split()[0])==istate:
                    outs.append(line.split())
                    istate = istate + 1
            except ValueError:
                continue
    if verb: print('Got spectrum from '+filename)
    return outs

def count_duplicates(spec,match_columns,round_columns=None): #adds last column to specrum counting duplicate states
    if round_columns is not None:
        for col_num in round_columns:
            spec[:,col_num] = [str(round_to_half(float(x))) for x in spec[:,col_num]]
    counts = times_so_far(list([list(tj) for tj in spec[:,match_columns]]))
    counts = np.array(counts)
    counts = counts.reshape(-1,1)
    return np.append(spec,counts,axis=1)

def cleanup():
    cleanup_cmd = " ".join(["rm","*.lcoef","timingdata.bigstick",\
        "timinginfo.bigstick","fort.*"])
    subprocess.Popen(cleanup_cmd,shell=True)

def remove_file(fn):
    if os.path.exists(fn):
        os.remove(fn)


#
#   SENSITIVITY ANALYSIS
#


def perturb_standard(pert_name,pert_vec):
    # pert_name has the  form usdb_p<idx>_<dv>
    remove_file(pert_name+".vec")
    remove_file(pert_name+".int")
    # PERTURB USDB VECTOR
    anamil_args = ["p","66",fn_int_vec_sorted]
    for pert in pert_vec:
        anamil_args.append(str(int(pert[0]))+" "+str(pert[1]))
    anamil_args.append("0 0")
    anamil_args.append(pert_name+".vec")
    anamil_args.append("x")
    arun = subprocess.Popen(anamil_cmd,stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    arun.stdin.write("\n".join(anamil_args).encode())
    arun_out = arun.communicate()[0]
    #arun.wait()
    arun.stdin.close()
    if verb: print("anamil args:\n"+"\n".join(anamil_args))

    # CONVERT TO INT
    standard_args = ["r","sd",pert_name,pert_name+".vec"]
    srun = subprocess.Popen(standard_cmd,stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    srun.stdin.write("\n".join(standard_args).encode())
    srun_out = srun.communicate()[0]
    #srun.wait()
    srun.stdin.close()
    if verb: print("standardorderint args:\n"+"\n".join(standard_args))

def perturb_milcom(pert_name,pert_vec):
    pert_mil_name = pert_name[:4]+'mil'+pert_name[4:]
    # pert_name has the  form usdb_m<idx>_<dv>
    remove_file(pert_name+".vec")
    remove_file(pert_mil_name+".vec")
    remove_file(pert_name+".int")
    # PERTURB MILCOMS VECTOR
    pert_mil_args = ["p","66",fn_milcoms_vector]
    for pert in pert_vec:
        pert_mil_args.append(str(int(pert[0]))+" "+str(pert[1]))
    pert_mil_args.append("0 0")
    pert_mil_args.append(pert_mil_name+".vec")
    pert_mil_args.append("x")
    pert_mil_run = subprocess.Popen(anamil_cmd,stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    pert_mil_run.stdin.write("\n".join(pert_mil_args).encode())
    pert_mil_run_out = pert_mil_run.communicate()[0]
    #pert_mil_run.wait()
    pert_mil_run.stdin.close()
    if verb: print("perturb milcom args:\n"+"\n".join(pert_mil_args))

    # TRANSFORM TO USDB BASIS
    trans_mil_args = ["t","66",fn_milcoms,pert_mil_name+'.vec','f'\
        ,pert_name+'.vec']
    trans_mil_args.append("x")
    trans_mil_run = subprocess.Popen(anamil_cmd,stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    trans_mil_run.stdin.write("\n".join(trans_mil_args).encode())
    trans_mil_run_out = trans_mil_run.communicate()[0]
    #trans_mil_run.wait()
    trans_mil_run.stdin.close()
    if verb: print("transform milcom args:\n"+"\n".join(trans_mil_args))

    # CONVERT TO INT
    standard_args = ["r","sd",pert_name,pert_name+".vec"]
    srun = subprocess.Popen(standard_cmd,stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    srun.stdin.write("\n".join(standard_args).encode())
    srun_out = srun.communicate()[0]
    #srun.wait()
    srun.stdin.close()
    if verb: print("standardorderint args:\n"+"\n".join(standard_args))

def compute_bigstick_spectrum(ZNpair,twoJz,filename_int,keep_wfn,fhlog=None):

    if fhlog is None:
        make_log = False
    else:
        make_log = True

    case_name = "Z"+str(int(ZNpair[0])).zfill(2)+"N"+str(int(ZNpair[1])).zfill(2)+'_2Jz'+str(int(twoJz))
    Zv = int(ZNpair[0]-8)
    Nv = int(ZNpair[1]-8)
    A = Zv+Nv+16

    if keep_wfn:
        b_opt = 'n'
    elif not keep_wfn:
        b_opt = 'ns'

    serial_condition = (twoJz>14 or any([Zv<2,Nv<2,Zv>10,Nv>10]))
    if serial_condition:
        if verb: print("RUNNING IN SERIAL")
        bigstick_cmd = bigstick_serial
        lanstring='ex'
        bargs = [b_opt,case_name,"sd"," ".join([str(Zv),str(Nv)]),str(twoJz),filename_int,\
            " ".join(["1","18",str(A),"0.3"]),"end",lanstring," ".join([str(nkeep),str(liter)])]
    else:
        if verb: print("RUNNING IN PARALLEL")
        bigstick_cmd = [bigstick_parallel]
        lanstring = "ld"
        bargs = [b_opt,case_name,"sd"," ".join([str(Zv),str(Nv)]),str(twoJz),filename_int,\
            " ".join(["1","18",str(A),"0.3"]),"end",lanstring," ".join([str(nkeep),str(liter)])]
    if verb: print(bargs)
    brun = subprocess.Popen(bigstick_cmd,stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    brun.stdin.write("\n".join(bargs).encode())
    bout = brun.communicate()[0]
    #brun.wait()
    bout = bout.decode('UTF-8')
    brun.stdin.close()
    with open(case_name+'.stdout','w') as boutfh:
        bout = bout.split('\n')
        boutfh.write("\n".join(bout))

    spectrum = np.array(getspectrum(case_name+".res")).astype(float)
    spectrum = count_duplicates(spectrum,match_columns=[3,4],round_columns=[3,4])
    dtype_list = [('#','int'),('E','float'),('Ex','float'),('J','float'),('T','float'),('n','int')]
    spectrum = np.array([tuple(x) for x in spectrum],dtype=dtype_list)
    if make_log:
        try:
            fhlog.write("test: "+str(spectrum[0])+"\n")
        except IndexError:
            print("BIGSTICK ERROR")
            print("\n".join(bargs))
            fhlog.write("\nBIGSTICK ERROR \n")
            fhlog.write("\n".join(bargs))
            exit()

        try:
            fhlog.write("test: "+str(spectrum[0])+"\n")
        except IndexError:
            print("BIGSTICK ERROR")
            print("\n".join(bargs))
            fhlog.write("\nBIGSTICK ERROR \n")
            fhlog.write("\n".join(bargs))
            exit()
    return spectrum

def match_states(current_data, spectrum, twoJz,fhlog=None):
    current_nstates = current_data.shape[0]

    if fhlog is None:
        make_log = False
    else:
        make_log = True

    #only keep states with 2*J = twoJz
    #spectrum = spectrum[ (2*spectrum['J']).astype(int) == twoJz]

    for i,current_data_state in enumerate(current_data):   #match states in current spectra
        for j,spectrum_state_tup in enumerate(spectrum):
            spectrum_state = np.array(list(spectrum_state_tup))
            current_data_state_2j2tn = np.array([float(x) for x in current_data_state[[3,4,5]]])
            current_data_state_2j2tn = current_data_state_2j2tn.astype(int)
            spectrum_state_2j2tn = np.array([float(x) for x in spectrum_state[[3,4,5]]])
            spectrum_state_2j2tn[[0,1]] = 2.0 * spectrum_state_2j2tn[[0,1]]
            spectrum_state_2j2tn = spectrum_state_2j2tn.astype(int)
            if (current_data_state_2j2tn == spectrum_state_2j2tn).all():
                if verb:
                    print("MATCHED")
                    print("exp:\t%i %i %i %i %i %i %f %f %f" % tuple(current_data[i]))
                    print("sm:\t"+" "+str(spectrum_state))
                if make_log:
                    fhlog.write("\nMATCHED \n")
                    fhlog.write("exp:\t%i %i %i %i %i %i %f %f %f \n" % tuple(current_data[i]))
                    fhlog.write("sm:\t"+str(spectrum_state)+"\n")

                current_data[i,-1] = spectrum_state[1]   # put energy in current_data

                spectrum = np.delete(spectrum,(j),axis=0)
                break
    return current_data, spectrum

def compute_energies(exp_data,filename_int,keep_wfn=False,make_log=False):
    error_states=[]
    if filename_int.endswith('.int'):
        filename_int = filename_int[:-4]
    if make_log:
        fhlog = open(filename_int+".log",'w')
    else:
        fhlog = None
    cases=[]  # cases defined by Z,N for individual bigstick run
    for pair in [list(p) for p in list(exp_data[:,1:3])]:
        if pair not in cases:
            cases.append(pair)
    cases = np.array(cases)
    ndata = len(exp_data)
    #Ebigstick = np.zeros(ndata) # absolute energies from bigstick corresponding to exp_data
    case_data_list = []
    for ncase,ZNpair in enumerate(cases):
        #case_name = "Z"+str(int(ZNpair[0]))+"N"+str(int(ZNpair[1]))
        if verb: print("COMPUTING CASE Z=%i N=%i" % tuple(ZNpair))
        if make_log: fhlog.write("\nCOMPUTING CASE Z=%i N=%i \n" % tuple(ZNpair))

        current_data = [] #exp data for case = ZNpair
        for i,entry in enumerate(exp_data):
            if not (entry[1:3]-ZNpair).any():
                current_data.append(entry)
        current_data = np.array(current_data)
        if make_log: fhlog.write("\nEXP DATA: \n")
        if verb:
            print("EXP DATA:")
            for current_line in current_data:
                print("%i %i %i %i %i %i %f %f" % tuple(current_line))
                if make_log: fhlog.write("%i %i %i %i %i %i %f %f \n" % tuple(current_line))

        # add last column to current_data for bigstick energies
        current_data = np.hstack((current_data,np.zeros((current_data.shape[0],1))))
        #current_data = np.hstack((current_data,-1*np.ones((current_data.shape[0],1))))

        debug = False     # BYPASS OBSERVABLE CALCULATION
        if not debug:
            twoJ_unique = np.unique(current_data[:,3]).astype(int)
            if verb: print(["2J unique:",twoJ_unique])
            for twoJz in twoJ_unique:
                if verb:
                    if any(current_data[current_data[:,3]==twoJz,-1]==0.0):
                    #if any(current_data[current_data[:,3]==twoJz,-1]<0.0):
                        print("COMPUTING WITH twoJz = %i" % twoJz)
                    else:
                        print("ALL STATES ACCOUNTED FOR WITH twoJz = %i" % twoJz)
                        continue
                spectrum = compute_bigstick_spectrum(ZNpair,twoJz,filename_int,keep_wfn,fhlog)

                if verb: print("BIGSTICK SPECTRUM:")
                if make_log: fhlog.write(" \nBIGSTICK SPECTRUM: \n")
                for spec_line in spectrum:
                    if verb: print(spec_line)
                    if make_log: fhlog.write(str(spec_line)+"\n")

                current_data, spectrum = match_states(current_data, spectrum, twoJz, fhlog)

            if not (current_data[:,-1]).all():
            #if (current_data[:,-1] < 0.0).any():
                print('MISSING STATE:')
                for current_line in current_data:
                    print("%i %i %i %i %i %i %f %f %f" % tuple(current_line))
                exit()

        case_data_list.append(current_data)

    output_data = np.vstack(case_data_list)
    #dtype_list = [('A','int'),('Z','int'),('N','int'),('2j','int'),('2t','int'),('n','int'),('Eexp','float'),('dEexp','float'),('Esm','float')]
    #output_data = np.array(output_data,dtype=dtype_list)
    fnout = filename_int +'_appended.dat'
    np.savetxt(fnout,output_data,fmt="\t".join(['%i']*6+['%f']*3))
    if make_log: fhlog.close()

def compute_chi_squared(filename_int,use_exp_errors):

    filename_in = filename_int +'_appended.dat'
    data = np.loadtxt(filename_in,delimiter="\t")
    print("DATA SHAPE =" + str(data.shape))
    if use_exp_errors == True:
        ndata = data.shape[0]
        exp_err_vec = np.sqrt( data[:,7]**2 + (theory_err*np.ones(ndata))**2 )
        residual = (data[:,-1] - data[:,6]) / exp_err_vec
    elif use_exp_errors == False:
        residual = (data[-1] - data[6])
    chisq = np.dot(residual,residual)

    if verb: print("MAX ABS DIFF = "+str(max([abs(ed) for ed in residual])))
    if verb: print("CHI SQUARED = "+str(chisq))
    return chisq



#
# EXPECTATION VALUES
#

def make_ops_vec(int_vec):
    nops = len(int_vec)
    for iop in range(nops):
        basis_vec = np.zeros(nops)
        basis_vec[iop] = 1
        #opvec = np.multiply(int_vec,basis_vec)
        outfn = interaction_name+'_o'+str(iop+1).zfill(2)+'.vec'
        np.savetxt(outfn,basis_vec,delimiter="\n",fmt="%f10",header=str(nops),comments='')

def make_ops_int(int_vec):
    nops = len(int_vec)
    for iop in range(nops):
        op_name = interaction_name+'_o'+str(iop+1).zfill(2)
        op_name_vec = op_name+'.vec'
        op_name_int = op_name+'.int'
        if os.path.isfile(op_name_int):
            rm_cmd = "rm "+op_name_vec
            subprocess.Popen(rm_cmd,shell=True)
        # CONVERT TO INT
        standard_args = ["r","sd",op_name,op_name_vec]
        srun = subprocess.Popen(standard_cmd,stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        srun.stdin.write("\n".join(standard_args).encode())
        srun_out = srun.communicate()[0]
        srun.stdin.close()
        #srun.wait()
        if verb: print("standardorderint args:\n"+"\n".join(standard_args))


def compute_exp_values(filename_int, make_log=False):
    #filename_int is operator .int file
    os.chdir(".")
    cases = glob.glob("Z*N*.wfn")  # this checks EVERY wavefunction present in the directory
    if filename_int.endswith('.int'):
        filename_int = filename_int[:-4]
    if make_log:
        fhlog = open(filename_int+".log",'w')
    else:
        fhlog = None
    op_name = filename_int[-3:]
    for icase,fnwfn in enumerate(cases):
        casename = fnwfn[:-4]
        if verb: print(casename)
        ZNpair = [int(casename[1:3]),int(casename[4:6])]
        casename_op = casename + "_" + op_name
        Zv = int(ZNpair[0]-8)
        Nv = int(ZNpair[1]-8)
        A = Zv+Nv+16
        if verb: print("COMPUTING CASE " + casename_op)
        if make_log: fhlog.write("\nCOMPUTING CASE " + casename_op +" \n ")
        ##
        lanstring = "ld"
        bigstick_cmd = bigstick_serial
        bargs = ["x",casename,casename_op,filename_int,\
               " ".join(["1","18",str(A),"0.3"]),"end"]
        #bargs = ["x",casename,casename_op,filename_int,\
        #        " ".join(["1","1","1","1"]),"end"]
        #if Zv<2 or Nv<2 or Zv>10 or Nv>10:
        #   bigstick_cmd = bigstick_name
        #   lanstring='ex'
        #   bargs = ["n",case_name,"sd"," ".join([str(Zv),str(Nv)]),str(A%2),filename_int,\
        #           " ".join(["1","18",str(A),"0.3"]),"end",lanstring," ".join([str(nkeep),str(liter)])]
        brun = subprocess.Popen(bigstick_cmd,stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        brun.stdin.write("\n".join(bargs).encode())
        bout = brun.communicate()[0]
        #brun.wait()
        brun.stdin.close()
        bout = str(bout)
        bout = bout.split('\n')
        with open(casename+'_x.stdout','w') as boutfh:
            boutfh.write("\n".join(bout))
        if make_log: fhlog.write("\n".join(bargs))
    if make_log: fhlog.close()


def gather_exp_values(exp_data,fhout,filename_int,make_log = False):
    op_num = int(filename_int[-2:])
    cases = glob.glob("Z*N*o"+str(op_num).zfill(2)+"*.res")
    cases = sorted(cases,key=operator.itemgetter(1,2))
    #print(cases)
    if make_log:
        fhlog = open(filename_int+"_gather.log",'w')
    else:
        fhlog = None
    #out_array = []
    match_list=[]
    for cn in cases:
        if make_log: fhlog.write("CASE "+cn+"\n")
        if verb: print("CASE "+cn)
        Zi = int(cn[1:3])
        Ni = int(cn[4:6])
        spectrum = np.array(getspectrum(cn))
        spectrum = spectrum.astype(np.float)
        spectrum = count_duplicates(spectrum,match_columns=[2,3],round_columns=None)
        for ed in exp_data:
            found = False
            if ("".join(str(ed)) in match_list):
                continue
            if ed[1]==Zi and ed[2]==Ni:
                for sd in spectrum:
                    #match T^2 = T(T+1)
                    exp_jtn = [ed[3],ed[4],ed[5]]   # exp J T n
                    spec_jtn = [sd[2],sd[3],sd[6]]  #spectrum J T n
                    exp_jtn[0] = float(exp_jtn[0])/2 #exp J
                    exp_jtn[1] = (float(exp_jtn[1])/2)*((float(exp_jtn[1])/2) + 1) #exp T^2
                    if verb:
                        print('exp_jtn',exp_jtn)
                        print('spec_jtn',spec_jtn)
                    if make_log: fhlog.write("CHECK jtn:\n"+str(exp_jtn)+"\n"+str(spec_jtn)+"\n")
                    if exp_jtn==spec_jtn:
                        if verb: print(sd)
                        # outline = Z   N   J   T   n   Esm    Eexp     dEexp   op#     <o>
                        outline = [str(Zi),str(Ni)]+["%i" % (x) for x in spec_jtn]\
                            +["%f" % x for x in [sd[1],ed[6],ed[7]]]+[str(op_num),str(sd[4])]
                        outstring = "\t".join(outline) + "\n"
                        fhout.write(outstring)
                        #out_array.append(outline)
                        if make_log: fhlog.write("MATCH \n"+str(ed)+"\n"+str(sd)+"\n")
                        if make_log: fhlog.write("OUT\n"+outstring)
                        if verb:
                            print('MATCH',ed,sd)
                            print('CASE',Zi,Ni)
                            print("OUT",outstring)
                        found = True
                        match_list.append("".join(str(ed)))
                        break
                    else:
                        if make_log: fhlog.write("NOPE \n"+str(ed)+"\n"+str(sd)+"\n")
                        #print('NO MATCH',ed,sd)
                if not found:
                    if make_log: fhlog.write("MISSING"+"\n"+str(ed)+"\n")

    #out_array =
    if make_log: fhlog.close()








