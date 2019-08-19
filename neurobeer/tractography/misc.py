""" misc.py

Module providing miscallaneous functionality

"""

import os.path as op
import numpy as np

def saveEig(dir_path, eigval_arr, eigvec_arr, verbose=0):
    """
    Function used to save eigenvalues and eigenvectors to binary file.

    INPUT:
        dir_path - directory path for storing files
        eigval_arr - array of eigenvalues
        eigvec_arr - matrix of eigenvectors

    OUTPUT:
        none
    """
    # Paths to save
    eigval_path = op.join(op.realpath(dir_path), "eigval.npy")
    eigvec_path = op.join(op.realpath(dir_path), "eigvec.npy")

    # Save file
    np.save(eigval_path, eigval_arr)
    np.save(eigvec_path, eigvec_arr)

    vprint("Saved eigenvalues & eigenvectors to %s" % op.realpath(dir_path),
            verbose)

    return

def vprint(txt, verbose):
    """
    Function used to print verbose statements

    INPUT:
        txt - message to return
        verbose - verbosity

    OUTPUT:
        none
    """
    if verbose != 0:
        print(txt)

    return
