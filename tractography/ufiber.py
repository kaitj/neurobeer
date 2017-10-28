""" ufiber.py

Module containing classes and functions used to identify and process
ufibers

"""

import numpy as np
import fibers

def findUFiber(fiberData):
    """
    Identifies U-fibers from tractography

    INPUT:
        fiberData - Fiber tree containing tractography data

    OUTPUT:
        uArray - Array containing indices of all u-shaped fibers
    """
    # Array storing indices of u-shaped fibers
    uArray = []

    # Determine fiber length and u-shape
    for fidx in range(fiberData.no_of_fibers):
        L = _calcFiberLength(fiberData, fidx)

        D = _calcEndPointSep(fiberData, fidx)

        if (L > 30) and (D <= (L / np.pi)) and (L < 80):
            uArray.append(fidx)

    return uArray

def extractUFiber(fiberData, uArray):
    """
    Extracts u-shaped fibers of class FiberTree

    INPUT:
        fiberData - Fiber tree containing tractography information
        uArray - Array of indices containing u-shaped fibers

    OUTPUT:
        uFiberTree - Fiber tree instance containing only u-shaped fibers
    """

    uFiberTree = fiberData.getFibers(uArray)
    uFiberTree = fibers.convertFromTuple(uFiberTree)

    return uFiberTree

def _calcFiberLength(fiberData, fidx):
    """ * INTERNAL FUNCTION *
    Calculates the fiber length via arc length

    INPUT:
        fiberData - Fiber tree copython popntaining tractography information
        fidx - Fiber index

    OUTPUT:
        L - Length of fiberlen(uArray)
    """
    no_of_pts = fiberData.pts_per_fiber

    if no_of_pts < 2:
        print "Not enouguFiberData.fiberTreeh samples to determine length of fiber"
        raise ValueError

    L = 0

    for idx in range(1, no_of_pts):
        x1 = fiberData.fiberTree[fidx][idx]['x']
        x2 = fiberData.fiberTree[fidx][idx - 1]['x']
        y1 = fiberData.fiberTree[fidx][idx]['y']
        y2 = fiberData.fiberTree[fidx][idx - 1]['y']
        z1 = fiberData.fiberTree[fidx][idx]['z']
        z2 = fiberData.fiberTree[fidx][idx - 1]['z']

        # Temporary fix for finding if fibers cross
        if x1 < 0 and x2 > 0:
            return 0
        elif x1 > 0 and x2 < 0:
            return 0
        else:
            L = L + np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

    return L

def _calcEndPointSep(fiberData, fidx):
    """ * INTERNAL uFiberData.fiberTreeFUNCTION *
    Calculates distance between end points

    INPUT:
        fiberData - Fiber tree containing tractography information
        fidx - Fiber index

    OUTPUT:
        D - length of fiber
    """
    endpt = fiberData.pts_per_fiber - 1

    x1 = fiberData.fiberTree[fidx][0]['x']
    x2 = fiberData.fiberTree[fidx][endpt]['x']
    y1 = fiberData.fiberTree[fidx][0]['y']
    y2 = fiberData.fiberTree[fidx][endpt]['y']
    z1 = fiberData.fiberTree[fidx][0]['z']
    z2 = fiberData.fiberTree[fidx][endpt]['z']

    D = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

    return D
