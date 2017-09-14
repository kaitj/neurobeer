""" misc.py

Module containing miscallaneous functions used in the library. These functions do not fit
in with other modules.

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def corr(data1, data2):
    """
    Perform correlation coefficient calculation using numpy library.
    Used for parallel, row-by-row correlation, returning only off-diagonal value.

    INPUT:
        data1 - First row of array data to be used in computation
        data2 - Second row of array data to be used in computation

    OUTPUT:
        val - Off-diagonal correlation value
   """

    val = np.corrcoef(data1, data2)
    val = val[0][1]

    return val

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
