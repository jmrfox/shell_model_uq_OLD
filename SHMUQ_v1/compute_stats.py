# compute stats for pkl file of B-values
import numpy as np
import pickle
import sys


save_hist = True

def print_stats(case_name,run_name,trans_idx):
    strengths = []
    for inter in data:
        strengths.append(data[inter][case_name][trans_idx]['B'])
    mu = np.mean(strengths)
    sigma = np.std(strengths)
    sig_unc = sigma/np.sqrt(2*len(strengths) - 2)
    #print("Transition index = %i \n mu = %f , 1-sigma = %f  +/-  %f " % (trans_idx,mu,sigma,sig_unc))
    print("mu = %f , 1-sigma = %f  +/-  %f " % (mu,sigma,sig_unc))
    if save_hist:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.ioff()
        plt.hist(strengths,bins=50)
        plt.yticks([])
        plt.savefig(run_name+'_trans'+str(trans_idx)+'.pdf')

if __name__=="__main__":
    fp = sys.argv[1]
    with open(fp,'rb') as fh:
        data = pickle.load(fh)

    fn = fp.split('/')[-1]
    case_name = fn.split('_')[0]
    run_name = fn.strip('.pkl')

    keys = data.keys()
    some_key = [x for x in keys][0]
    #n_trans = len(data['usdb_rand0000'][case_name])
    n_trans = len(data[some_key][case_name])
    for trans_idx in range(n_trans):
        op_type = data[some_key][case_name][trans_idx]['type']
        ji = data[some_key][case_name][trans_idx]['2Ji']//2
        ni = data[some_key][case_name][trans_idx]['ni']
        jf = data[some_key][case_name][trans_idx]['2Jf']//2
        nf = data[some_key][case_name][trans_idx]['nf']
        print("%s : %i(%i)  --> %i(%i) " % (op_type,ji,ni,jf,nf))
        print_stats(case_name,run_name,trans_idx)


