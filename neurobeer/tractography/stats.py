""" stats.py

Module containing classes and functions used to compute quantitative
tract-based statistics.

"""

import os, csv
import numpy as np

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
        clusterAvg = np.nanmean(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                            scalarType)[:, :], axis=0)
        avg = np.nanmean(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                      scalarType)[:, :])
    else:
        clusterAvg = np.nanmean(fiberTree.getScalars(idxes, scalarType)[:, :],
                             axis=0)
        avg = np.nanmean(fiberTree.getScalars(idxes, scalarType)[:, :])

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
        clusterSdev = np.nanstd(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                             scalarType)[:, :], axis=0)
        stdev = np.nanstd(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                       scalarType)[:, :])
    else:
        clusterSdev = np.nanstd(fiberTree.getScalars(idxes, scalarType)[:, :],
                             axis=0)
        stdev = np.nanstd(fiberTree.getScalars(idxes, scalarType)[:, :])
    return clusterSdev, stdev

def calcGeoStats(LArray):
    """
    Calculates the mean and standard deviation fiber length for an identified
    group of fibers

    INPUT:
        LArray - array of fiber lengths

    OUTPUT:
        LMean - mean fiber length
        LSD - standard deviation of fiber length
        fiberCount - number of fibers
    """

    LMean = np.nanmean(LArray)
    LSD = np.nanstd(LArray)
    fiberCount= len(LArray)

    return LMean, LSD, fiberCount

def writeGeoCSV(clusterLabel, LMean, LStd, fiberCount, dirpath=None):
    """
    Writes the length and distance of each cluster for an identified group of
    fibers

    INPUT:
        clusterLabel - label for cluster
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
        header = ['Cluster ID', 'Length Mean', 'Length S.D.', 'Fiber Count']
        writer = csv.DictWriter(f, delimiter=',', lineterminator='\n',
                                fieldnames=header)
        if not infoExists:
            writer.writeheader()
        writer = csv.writer(f)
        writer.writerow([clusterLabel, LMean, LStd, fiberCount])

    f.close()

def writeCSV(clusterLabel, fiberTree, scalarType, idxes=None, dirpath=None):
    """
    Writes stats to a csv file. Outputs one file per scalar type
    Each row of the csv is the average value at each sampled point

    INPUT:
        clusterLabel - label for cluster
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
        scalarArray = \
            np.nanmean(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                    scalarType)[:, :], axis=0)
    else:
        scalarArray = np.nanmean(fiberTree.getScalars(idxes, scalarType)[:, :],
                      axis=0)

    fileName = scalarType.split('/', -1)[-1]
    filePath = dirpath + '/stats/' + fileName + '.csv'

    fileExists = os.path.isfile(filePath)

    with open(filePath, 'a') as f:
        headers = ['Cluster ID']
        for idx in range(fiberTree.pts_per_fiber):
            headers.append(idx)
        writer = csv.DictWriter(f, delimiter=',', lineterminator='\n',
                                fieldnames=headers)
        if not fileExists:
            writer.writeheader()

        writer = csv.writer(f)
        rowInfo = [clusterLabel]
        for idx in range(fiberTree.pts_per_fiber):
            rowInfo.append(scalarArray[idx])
        writer.writerow(rowInfo)

    f.close()
