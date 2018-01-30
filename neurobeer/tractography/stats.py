""" stats.py

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

def calcGeoStats(LArray):
    """
    Calculates the mean and standard deviation fiber length for an identified group of fibers

    INPUT:
        LArray - array of fiber lengths

    OUTPUT:
        LMean - mean fiber length
        LSD - standard deviation of fiber length
        fiberCount - number of fibers
    """

    LMean = np.mean(LArray)
    LSD = np.std(LArray)
    fiberCount= len(LArray)

    return LMean, LSD, fiberCount

def writeGeoCSV(LMean, LStd, fiberCount, dirpath=None):
    """
    Writes the length and distance of each cluster for an identified group of fibers

    INPUT:
        LMean - mean length of cluster
        LStd - standard deviation of cluster
        fiberCount - number of fibers in cluster
        dirpath - directory to store CSV file; default None

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

    infoName = 'clusterInfo'
    infoPath = statspath + infoName + '.csv'
    infoExists = os.path.isfile(infoPath)

    with open(infoPath, 'a') as f:
        header = ['Length Mean', 'Length S.D.', 'Fiber Count']
        writer = csv.DictWriter(f, delimiter=',', lineterminator='\n', fieldnames=header)
        if not infoExists:
            writer.writeheader()
        writer = csv.writer(f)
        writer.writerow([LMean, LStd, fiberCount])

    f.close()

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
