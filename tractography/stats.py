""" stats.py

Module containing classes and functions used to compute quantitative
tract-based statistics.

"""

import os
import numpy as np
import matplotlib.pyplot as plt

def _mean(fiberTree, scalarType, idxes=None):
    """
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
    """
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

def plotStats(fiberTree, scalarType, idxes=None, dirpath=None):
    """
    Plots the calculated tract-based statistics for each fiber bundle

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        scalarType - type of quantitative data to plot
        idxes - indices to extract info from; defaults None (returns data for all fibers)
        dirpath - location to st#a = np.array([1, 1, 1], [2, 2, 2])ore plots; defaults None

    OUTPUT:
        none
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
    ytitle = scalarType.split('_', -1)[-1]

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
