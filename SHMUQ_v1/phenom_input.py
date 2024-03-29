#
#   Generates input for phenomenological calculation in bigstick
#   J Fox 2019
#

nl = "\n"
sp = " "

def make_input_dict():
    input_dict = {}
    input_dict['run_opt'] = None
    input_dict['run_name'] = None
    input_dict['sps_file'] = None
    input_dict['Zv'] = None
    input_dict['Nv'] = None
    input_dict['twoJz'] = None
    input_dict['frag_opt'] = None
    input_dict['lanc_opt'] = None
    input_dict['nkeep'] = None
    input_dict['lanc_iters'] = None
#     for i in range(n_ints)
#         input_dict['int_file_'+i] = None
#         input_dict['int_scale_'+i] = None
    return input_dict

def write_out_file(fn,d,n_ints):
    with open(fn,'w') as fh:
        fh.write(d['run_opt']+nl)
        fh.write(d['run_name']+nl)
        fh.write(d['sps_file']+nl)
        fh.write(d['Zv'] + sp +d['Nv']+nl)
        fh.write(d['twoJz']+nl)
        if d['frag_opt'] is not None:
            fh.write(d['frag_opt']+nl)
        for i in range(n_ints):
            fh.write(d['int_file_'+str(i)]+nl)
            fh.write(d['int_scale_'+str(i)]+nl)
        fh.write('end'+nl)
        fh.write(d['lanc_opt']+nl)
        fh.write(d['nkeep'] + sp + d['lanc_iters']+nl)


if __name__=="__main__":

    input_dict = make_input_dict()

#  MPI?
    x = input("Do you want MPI? (y/n)"+nl)
    if 'y' in x:
        input_dict['frag_opt'] = '0'

# option
    opt_list = ['n','d','dp']
    x = input("Enter bigstick option"+nl+str(opt_list)+nl)
    while x not in opt_list:
        x = input("Option not available. Try again.")
    input_dict['run_opt'] = x

# name
    x = input("Enter name for run"+nl)
    input_dict['run_name'] = x

# sps
    x = input("Enter name for .sps file"+nl)
    input_dict['sps_file'] = x

# Zv Nv
    x = input("Enter Zv Nv"+nl).split()
    input_dict['Zv'] = x[0]
    input_dict['Nv'] = x[1]

# twoJz
    x = input("Enter 2*Jz"+nl)
    input_dict['twoJz'] = x

# number of ints
    n_ints = int(input("Enter number of interactions"+nl))
    for i in range(n_ints):
        input_dict['int_file_'+str(i)] = input("Enter name of .int"+nl)
        input_dict['int_scale_'+str(i)] = \
                input("Enter scaling for .int   (i.e. 1 Acore+2 A 0.33) "+nl)

# lanczos
    x = input("Enter # of states"+nl)
    input_dict["lanc_opt"] = "ld"
    input_dict['nkeep'] = x
    input_dict["lanc_iters"] = str(20*int(x))

# out
    #fn_out = input("Enter input filename:"+nl)
    fn_out = input_dict['run_name']+'.input'
    write_out_file(fn_out,input_dict,n_ints)





