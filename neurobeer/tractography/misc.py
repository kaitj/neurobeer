""" misc.py

Module containing miscallaneous functions used in the library. These functions do not fit
in with other modules.

"""

import numpy as np

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
    NOTE: Matplotlib not memory-friendly, save matrix to textfile

    INPUT:
        dirpath - Directory path for storing matrix images
        matrix - Similarity matrix to be saved
        matrixType - Type of similarity matrix (ie. weighted, FA, T1)

    OUTPUT:
        none
    """

    fname = (dirpath + matrixType + '_Similarity.txt')

    np.savetxt(fname, matrix)

    # OLD CODE FOR PLOTTING FIGURE
    # f = plt.figure(figsize=(10, 10))
    # matrix = plt.imshow(matrix, cmap='viridis', interpolation='none')
    # plt.title((matrixType + ' Similarity'), fontsize=16)
    #
    # ax = plt.gca()
    # ax.tick_params(axis='both', labelsize=14)
    # for axis in [ax.xaxis, ax.yaxis]:
    #     axis.set_major_locator(ticker.MaxNLocator(integer=True))
    # div = make_axes_locatable(ax)
    # cax = div.append_axes('right', size='5%', pad=0.25)
    # cax.tick_params(labelsize=14)
    # plt.colorbar(matrix, cax)
    # plt.tight_layout()
    #
    # plt.savefig(dirpath + '/' + matrixType + '_Similarity.png')
    # plt.close(f)
    #
    # del matrix
