""" misc.py

Module containing miscallaneous functions used in the library. These functions do not fit
in with other modules.

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def saveEig(dirpath, eigvalArray, eigvecArray):
    """
    Function used to save eigenvalues and eigenvectors to binary file):

    INPUT:
        dirpath - Directory path for storing eigenvalues & eigenvectors
        eigvalArray - Array of eigenvalues to be saved
        eigvecArray - Matrix of eigenvectos to be saved

    OUTPUT:
        none
    """
    # Paths for files to be saved
    eigvalPath = dirpath + '/eigval.npy'
    eigvecPath = dirpath + '/eigvec.npy'

    # Save to file
    np.save(eigvalPath, eigvalArray)
    np.save(eigvecPath, eigvecArray)

def saveMatrix(dirpath, matrix, matrixType):
    """
    Function used to save similarity matrices.

    INPUT:
        dirpath - Directory path for storing matrix images
        matrix - Similarity matrix to be saved
        matrixType - Type of similarity matrix (ie. weighted, FA, T1)

    OUTPUT:
        none
    """

    f = plt.figure(figsize=(10, 10))
    im = plt.imshow(matrix, cmap='viridis')
    plt.title((matrixType + ' Similarity'), fontsize=16)

    ax = plt.gca()
    ax.tick_params(axis='both', labelsize=14)
    div = make_axes_locatable(ax)
    cax = div.append_axes('right', size='5%', pad=0.25)
    cax.tick_params(labelsize=14)
    plt.colorbar(im, cax)

    plt.savefig(dirpath + '/' + matrixType + '_Similarity.png')
    plt.close(f)
