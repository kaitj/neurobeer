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
    LArray = []
    DArray = []

    # Determine fiber length and u-shape
    for fidx in range(fiberData.no_of_fibers):
        L = _calcFiberLength(fiberData, fidx)

        D = _calcEndPointSep(fiberData, fidx)

        # Temporary max length constraint
        if (L > 30) and (D <= (L / np.pi)) and (L < 90):
            uArray.append(fidx)
            LArray.append(L)
            DArray.append(D)

    return uArray, LArray, DArray

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

def uFiberStats(LArray, DArray, fidxes, no_of_fibers=20):
    """
    Calculates the mean and standard deviation for fiber length and distance between
    end points for a group of fibers.

    INPUT:
        LArray - array of fiber lengths
        DArray - array of distances between end points for fibers
        fidxes - array of indices of fibers to calculate
        no_of_fibers - number of fibers in tractography; defaults 20

    OUTPUT:
        LMean - mean fiber length
        LSD - standard deviation of fiber length
        DMean - mean distance between end points
        DSD - standard deviation between end points
    """

    Ltemp = []
    Dtemp = []

    for fidx in range(no_of_fibers):
        if fidx in fidxes:
            Ltemp.append(LArray[fidx])
            Dtemp.append(DArray[fidx])

    LMean = np.mean(Ltemp)
    LSD = np.std(Ltemp)
    DMean = np.mean(Dtemp)
    DSD = np.std(Dtemp)

    return LMean, LSD, DMean, DSD

def _calcFiberLength(fiberData, fidx):
    """ * INTERNAL FUNCTION *
    Calculates the fiber length via arc length

    INPUT:
        fiberData - Fiber tree copython popntaining tractography information
        fidx - Fiber index

    OUTPUT:
        L - Length of fiber
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
