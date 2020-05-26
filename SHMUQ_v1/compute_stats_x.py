# compute styats for pkl file of expectation-values
import numpy as np
import pickle
import sys

fn = sys.argv[1]
with open(fn,'rb') as fh:
    data = pickle.load(fh)
    mu = np.mean(data)
    sigma = np.std(data)
    sig_unc = sigma/np.sqrt(2*len(data) - 2)
    print("<O>: mu = %f , 1-sigma = %f  +/-  %f " % (mu,sigma,sig_unc))

save_hist = True
if save_hist:
    import matplotlib.pyplot as plt
    plt.hist(data,bins=50)
    #plt.yticks([])
    mu= np.mean(data)
    sigma = np.std(data)
    plt.title("mean = %f , std. dev. = %f" % (mu,sigma))
    plt.ylabel('Counts')
    plt.xlabel(rf"$\langle l \cdot s \rangle$")
    plt.savefig('new_figure.pdf')
