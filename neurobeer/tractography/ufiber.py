""" ufiber.py

Module containing classes and functions used to identify and process
ufibers

"""

import os, csv
import numpy as np
from . import fibers

def findUFiber(fiberData):
    """
    Identifies U-fibers from tractography

    INPUT:
        fiberData - fiber tree containing tractography data

    OUTPUT:
        uArray - array containing indices of all u-shaped fibers
        LArray - array containing lengths of all u-shaped fibers
        DArray - array containing end point seperation distance
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
        if (L > 20) and (D <= (L / np.pi)) and (L < 80):
            uArray.append(fidx)
            LArray.append(L)
            DArray.append(D)

    return uArray, LArray, DArray

def _mean(fiberTree, scalarType, idxes=None):
    """ *INTERNAL FUNCTION*
    Finds the average of all fibers in bundle at specific sample points

    INPUT:
        fiberTree - tree containing spatial and quantitative information of
                    fibers
        scalarType - type of quantitative data to find average of
        idxes - indices to extract info from; defaults None (returns data for
                all fibers)

    OUTPUT:
        clusterAvg - calculated tract-based average
        avg - calculated average of data for group
    """

    if idxes is None:
        clusterAvg = np.mean(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                             scalarType)[:, :], axis=0)
        avg = np.mean(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                      scalarType)[:, :])
    else:
        clusterAvg = np.mean(fiberTree.getScalars(idxes, scalarType)[:, :],
                             axis=0)
        avg = np.mean(fiberTree.getScalars(idxes, scalarType)[:, :])

    return clusterAvg, avg

def _stddev(fiberTree, scalarType, idxes=None):
    """ *INTERNAL FUNCTION*
    Finds the standard deviation of all fibers in bundle at specific sample
    points

    INPUT:
        fiberTree - tree containing spatial and quantitative information of
                    fibers
        scalarType - type of quantitative data to find standard deviation of
        idxes - indices to extract info from; defaults None (returns data for
                     all fibers)

    OUTPUT:
        clusterSdev - calculated tract-based standard deviation
        stdev - standard deviation of fiber group
    """
    if idxes is None:
        clusterSdev = np.std(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                             scalarType)[:, :], axis=0)
        stdev = np.std(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                       scalarType)[:, :])
    else:
        clusterSdev = np.std(fiberTree.getScalars(idxes, scalarType)[:, :],
                             axis=0)
        stdev = np.std(fiberTree.getScalars(idxes, scalarType)[:, :])

    return clusterSdev, stdev

def extractUFiber(fiberData, uArray):
    """
    Extracts u-shaped fibers of class FiberTree

    INPUT:
        fiberData - fiber tree containing tractography information
        uArray - array of indices containing u-shaped fibers

    OUTPUT:
        uFiberTree - fiber tree instance containing only u-shaped fibers
    """

    uFiberTree = fiberData.getFibers(uArray)
    uFiberTree = fibers.convertFromTuple(uFiberTree)

    return uFiberTree

def writeCSV(clusterLabel, LMean, LStd, DMean, DStd, fiberCount, dirpath=None):
    """
    Writes the length and distance of each cluster for a group of fibers.

    INPUT:
        clusterLabel - label for cluster
        LMean - mean of length for a cluster
        LStd - standard deviation of length of a cluster
        DMean - mean of distance between end points for a cluster
        DStd - standard deviation of distance between end points for a cluster

    OUTPUT:
        none
    """

    if dirpath is None:
        dirpath = os.getcwd()
    else:
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

    statspath = dirpath + '/stats/'
    if not os.path.exists(statspath):
        os.makedirs(statspath)

    fileName = 'clusterInfo'
    filePath = dirpath + '/stats/' + fileName + '.csv'
    fileExists = os.path.isfile(filePath)

    with open(filePath, 'a') as f:
        header = ['Cluster ID', 'Length Mean', 'Length S.D.', 'Distance Mean',
                  'Distance S.D.', 'Fiber Count']
        writer = csv.DictWriter(f, delimiter=',', lineterminator='\n',
                                fieldnames=header)
        if not fileExists:
            writer.writeheader()
        writer = csv.writer(f)
        writer.writerow([clusterLabel, LMean, LStd, DMean, DStd, fiberCount])

    f.close()

def uFiberStats(LArray, DArray, fidxes):
    """
    Calculates the mean and standard deviation for fiber length and distance
    between end points for a group of fibers.

    INPUT:
        LArray - array of fiber lengths
        DArray - array of distances between end points for fibers
        fidxes - array of indices of fibers to calculate

    OUTPUT:
        LMean - mean fiber length
        LSD - standard deviation of fiber length
        DMean - mean distance between end points
        DSD - standard deviation between end points
    """

    Ltemp = []
    Dtemp = []

    for fidx in range(len(LArray)):
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
        fiberData - fiber tree containing tractography information
        fidx - fiber index

    OUTPUT:
        L - fength of fiber
    """
    no_of_pts = fiberData.pts_per_fiber

    if no_of_pts < 2:
        raise ValueError("Not enough samples to determine length of fiber")

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
    """ * INTERNAL FUNCTION *
    Calculates distance between end points

    INPUT:
        fiberData - fiber tree containing tractography information
        fidx - fiber index

    OUTPUT:
        D - distance between end points
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
