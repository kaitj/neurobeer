""" ufiber.py

Module containing classes and functions used to identify and process
ufibers

"""

import os, csv
import numpy as np
import fibers

def findUfiber(fiberData):
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
            DArray.append(D)""" stats.py

Module containing classes and functions used to compute quantitative
tract-based statistics.

"""

import os, csv
import numpy as np
import matplotlib.pyplot as plt

def _mean(fiberTree, scalarType, idxes=None):
    """ *INTERNAL FUNCTION*
    Finds the average of all fibers in bundle at specific sample points

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        scalarType - type of quantitative data to find average of
        idxes - indices to extract info from; defaults None (returns data for all fibers)

    OUTPUT:
        clusterAvg - calculated tract-based average
        avg - calculated average of data for group
    """

    if idxes is None:
        clusterAvg = np.mean(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                            scalarType)[:, :], axis=0)
        avg = np.mean(fiberTree.getScalars(range(fiberTree.no_of_fibers), scalarType)[:, :])
    else:
        clusterAvg = np.mean(fiberTree.getScalars(idxes, scalarType)[:, :], axis=0)
        avg = np.mean(fiberTree.getScalars(idxes, scalarType)[:, :])

    return clusterAvg, avg

def _stddev(fiberTree, scalarType, idxes=None):
    """ *INTERNAL FUNCTION*
    Finds the standard deviation of all fibers in bundle at specific sample points

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        scalarType - type of quantitative data to find standard deviation of
        idxes - indices to extract info from; defaults None (returns data for all fibers)

    OUTPUT:
        clusterSdev - calculated tract-based standard deviation
        stdev - standard deviation of fiber group
    """
    if idxes is None:
        clusterSdev = np.std(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                            scalarType)[:, :], axis=0)
        stdev = np.std(fiberTree.getScalars(range(fiberTree.no_of_fibers), scalarType)[:, :])
    else:
        clusterSdev = np.std(fiberTree.getScalars(idxes, scalarType)[:, :], axis=0)
        stdev = np.std(fiberTree.getScalars(idxes, scalarType)[:, :])
    return clusterSdev, stdev

def writeCSV(fiberTree, scalarType, idxes=None, dirpath=None):
    """
    Writes stats to a csv file. Outputs one file per scalar type
    Each row of the csv is the average value at each sampled point

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        scalarType - type of quantitative data to plot
        idxes - indices to extract info from; defaults None (returns data for all fibers)
        dirpath - location to store plots; defaults None

    OUTPUT:
        none
    """

    if dirpath is None:
        dirpath = os.getcwd()
    else:
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

    if idxes is None:
        scalarArray = np.mean(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                            scalarType)[:, :], axis=0)
    else:
        scalarArray = np.mean(fiberTree.getScalars(idxes, scalarType)[:, :], axis=0)

    fileName = scalarType.split('/', -1)[-1]
    filePath = dirpath + '/stats/' + fileName + '.csv'

    fileExists = os.path.isfile(filePath)

    with open(filePath, 'a') as f:
        headers = range(1, fiberTree.pts_per_fiber+1)
        writer = csv.DictWriter(f, delimiter=',', lineterminator='\n', fieldnames=headers)
        if not fileExists:
            writer.writeheader()

        writer = csv.writer(f)
        writer.writerow(scalarArray)

    f.close()

def plotStats(fiberTree, scalarType, idxes=None, dirpath=None):
    """
    Plots the calculated tract-based statistics for each fiber bundle

    INPUT:np.mean(fiberTree.getScalars(range(fiberTree.no_of_fibers), scalarType)[:, :])
        fiberTree - tree containing spatial and quantitative information of fibers
        scalarType - type of quantitative data to plot
        idxes - indices to extract info from; defaults None (returns data for all fibers)
        dirpath - location to store plots; defaults None

    OUTPUT:
        nonenp.mean(fiberTree.getScalars(range(fiberTree.no_of_fibers), scalarType)[:, :])
    """

    if dirpath is None:
        dirpath = os.getcwd()
    else:
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

    # Statistical calculations for plot
    x = range(fiberTree.pts_per_fiber)
    yavg, avg = _mean(fiberTree, scalarType, idxes)
    ystd, stdev = _stddev(fiberTree, scalarType, idxes)

    # Plot of stats
    f = plt.figure(figsize=(10, 10))
    plt.grid()

    avgline = plt.plot(x, yavg, 'b', linewidth=3, label='Mean')
    avgdot = plt.plot(x, yavg, '.r', markersize=15)
    stdupper = plt.plot(x, yavg+ystd, 'k', linewidth=2, alpha=0.7)
    stdlower = plt.plot(x, yavg-ystd, 'k', linewidth=2, alpha=0.7, label='Std. Dev')
    plt.fill_between(x, yavg-ystd, yavg+ystd, facecolor='grey', alpha=0.7)
    plt.xlim(min(x), max(x))

    # Plot labels
    fileName = scalarType.split('/', -1)[-1]
    title = fileName + ' (%.2f +/- %.2f)' % (avg, stdev)
    ytitle = fileName

    plt.title(title, fontsize=16)
    plt.ylabel(ytitle.upper(), fontsize=14)
    plt.xlabel('Fiber Samples', fontsize=14)
    plt.legend(loc='upper left', fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=14,
                              length=10, pad=10, width=1, direction='out',
                              top='off', right='off')

    # Save figure
    plt.savefig(dirpath + '/' + fileName+ '_stats.png')
    plt.close(f)

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

def writeCSV(LMean, LStd, DMean, DStd, dirpath=None):
    """
    Writes the length and distance of each cluster for a group of fibers.

    INPUT:
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

    lengthName = 'clusterLengths'
    lengthPath = dirpath + '/stats/' + lengthName + '.csv'
    lengthExists = os.path.isfile(lengthPath)
    distanceName = 'clusterDistances'
    distancePath = dirpath + '/stats/' + distanceName + '.csv'
    distanceExists = os.path.isfile(distancePath)

    with open(lengthPath, 'a') as f1:
        header = ['Length Mean', 'Length S.D.']
        writer = csv.DictWriter(f1, delimiter=',', lineterminator='\n', fieldnames=header)
        if not lengthExists:
            writer.writeheader()
        writer = csv.writer(f1)
        writer.writerow([LMean, LStd])

    with open(distancePath, 'a') as f2:
        header = ['Distance Mean', 'Distance S.D.']
        writer = csv.DictWriter(f2, delimiter=',', lineterminator='\n', fieldnames=header)
        if not distanceExists:
            writer.writeheader()
        writer = csv.writer(f2)
        writer.writerow([DMean, DStd])

    f1.close()
    f2.close()

def uFiberStats(LArray, DArray, fidxes):
    """
    Calculates the mean and standard deviation for fiber length and distance between
    end points for a group of fibers.

    INPUT:
        LArray - array of fiber lengths
        DArray - array of distances between end points for fibers
        fidxes - array of indices of fibers to calculate
        no_of_fibers - number of fibers in tractography

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
