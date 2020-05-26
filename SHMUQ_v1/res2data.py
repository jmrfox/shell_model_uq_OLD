#
#   Convert a .res file from normal run to a data file matching the format of usdb_data.csv
#   Fox 4/20
#

import numpy as np
#import sys
from sa_mod import getspectrum,count_duplicates
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Convert a .res file from normal run to a data file matching the format of usdb_data.csv")

parser.add_argument('protons',metavar='Z',type=int,nargs=1,help="proton number")
parser.add_argument('neutrons',metavar='N',type=int,nargs=1,help="neutron number")

parser.add_argument('input',metavar='input',type=str,nargs=1,help=".res file to convert")
parser.add_argument('output',metavar='output',type=str,nargs=1,help=".txt file to print")

args = parser.parse_args()

fn_in = args.input[0]
s = getspectrum(fn_in)
s = np.array(s)
s = count_duplicates(s,match_columns=[3,4],round_columns=[3,4])
s =s.astype(float)

s[:,3] = 2*s[:,3]
s[:,4] = 2*s[:,4]

Z = args.protons[0]
N = args.neutrons[0]
A = Z+N

df = pd.DataFrame(s,columns=['i','energy','ex','2j','2t','n'])
df.insert(0,'N',N)
df.insert(0,'Z',Z)
df.insert(0,'A',A)
df.insert(0,'error',0.)

df.to_csv(args.output[0],columns=['A','Z','N','2t','2j','n','energy','error'],index=False)


print("Output written to %s" % args.output[0])
