""" misc.py

Module containing functions used in the tool. These functions do not fit
in other modules.

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def corr(data1, data2):
    """ Perform correlation coefficient calculation using numpy library.
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

def saveMatrix(dirpath, matrix, matrixType, colormap='viridis'):
    """ Function used to save matrices.
    """

    f = plt.figure(figsize=(10, 10))
    im = plt.imshow(matrix, cmap='viridis')
    plt.title((matrixType + 'Similarity'), fontsize=16)

    ax = plt.gca()
    ax.tick_params(axis='both', labelsize=14)
    div = make_axes_locatable(ax)
    cax = div.append_axes('right', size='5%', pad=0.25)
    cax.tick_params(labelsize=14)
    plt.colorbar(im, cax)

    plt.savefig(dirpath + '/' + matrixType.lower() + 'Similarity.png')
