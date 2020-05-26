#!/home/jfox/anaconda3/bin/python
#
#   Compute expectation values for given data
#	Jordan Fox 8/19 SDSU

from sa_mod import *
import subprocess
import argparse

parser = argparse.ArgumentParser(description="Computes expectation values from an input data set")
parser.add_argument('data_file',type=str,nargs=1,help="data file")
parser.add_argument('interaction_file',type=str,nargs=1,help="interaction file")
parser.add_argument('output',type=str,nargs=1,help="output file name for exp values")

args = parser.parse_args()

#fn_int = "usdb_sorted"
fn_data = args.data_file[0]
fn_int = args.interaction_file[0]

data, int_vec = get_data(fn_data)

#compute_energies(data,fn_int,keep_wfn=True)

make_ops_vec(int_vec)
make_ops_int(int_vec)

op_files = glob.glob(interaction_name + "_o??.int")
for op_fn in op_files:
    compute_exp_values(op_fn)

#fn_out = "exp_values_ne21data.dat"
fn_out = args.output[0]
with open(fn_out,'w') as fh_out:
    for op_fn in op_files:
        op_name = op_fn[:-4]
        gather_exp_values(data,fh_out,op_name)
