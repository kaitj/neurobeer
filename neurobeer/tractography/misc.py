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
    eigval_path = op.join(op.realpath(dir_path), "eigval.npz")
    eigvec_path = op.join(op.realpath(dir_path), "eigvec.npz")

    # Save file
    np.savez_compressed(eigval_path, eigval_arr)
    np.savez_compressed(eigvec_path, eigvec_arr)

    vprint("Saved eigenvalues & eigenvectors to %s" % op.realpath(dir_path),
            verbose)

    return

def vprint(txt, verbose, debug=False):
    """
    Function used to print verbose statements

    INPUT:
        txt - message to return
        verbose - verbosity
        debug - flag to indicate if message should be a debug message

    OUTPUT:
        none
    """
    if verbose == 1:
        print(txt)
    elif verbose == 3 and debug is True:
        print("DEBUG: %s" % txt)

    return
