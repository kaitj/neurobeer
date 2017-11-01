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
        avg - calculated average of data
    """

    if idxes is None:
        avg = np.mean(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                            scalarType)[:, :], axis=0)
    else:
        avg = np.mean(fiberTree.getScalars(idxes, scalarType)[:, :], axis=0)

    return avg

def _stddev(fiberTree, scalarType, idxes=None):
    """
    Finds the standard deviation of all fibers in bundle at specific sample points

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        scalarType - type of quantitative data to find standard deviation of
        idxes - indices to extract info from; defaults None (returns data for all fibers)

    OUTPUT:
        sdev - calculated standard deviation of data
    """
    if idxes is None:
        sdev = np.std(fiberTree.getScalars(range(fiberTree.no_of_fibers),
                            scalarType)[:, :], axis=0)
    else:
        sdev = np.std(fiberTree.getScalars(idxes, scalarType)[:, :], axis=0)
    return sdev

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

    # Info for plot labels
    title = scalarType.split('/', -1)[-1]
    title = title + '(%.2f +/- %.2f)' % np.mean(_mean, axis=1), np.mean(_stddev, axis=1)
    ytitle = scalarType.split('_', -1)[-1]

    # Statistical calculations for plot
    x = range(fiberTree.pts_per_fiber)
    yavg = _mean(fiberTree, scalarType, idxes)
    ystd = _stddev(fiberTree, scalarType, idxes)

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
    plt.title(title, fontsize=16)
    plt.ylabel(ytitle.upper(), fontsize=14)
    plt.xlabel('Fiber Samples', fontsize=14)
    plt.legend(loc='upper left', fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=14,
                              length=10, pad=10, width=1, direction='out',
                              top='off', right='off')

    # Save figure
    plt.savefig(dirpath + '/' + title + '_stats.png')
    plt.close(f)
