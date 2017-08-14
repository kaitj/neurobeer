""" misc.py

Module containing functions used in the tool. These functions do not fit
in other modules.

"""

import numpy as np

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
