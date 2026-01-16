import numpy as np

def get_ndim(ibzwfs):

    return ibzwfs.nbands

def get_n_occ(ibzwfs):

    f_n = ibzwfs.get_eigs_and_occs(kpt)[1]
    indices = f_n > 1.0e-5  #the third column contains the occupation numbers
    n_occ = sum(indices)
    return n_occ