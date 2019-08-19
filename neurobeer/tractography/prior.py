""" prior.py

Module containing classes and functions pertaining to use of previously clustered data.

"""

import os
import numpy as np
from . import fibers, tractio, misc
from vtk.util import numpy_support

def load(priorVTKPath, templateFlag=False, verbose=0):
    """
    Class used to load .vtk prior file.

    INPUT:
        priorVTKPath - absolute path to VTK file containing prior information
                       to be used
        templateflag - flag to set for subsetting; defaults false
        verbose - verbosity of function; defaults 0

    OUTPUT:
        priorTree - returns prior information stored in a tree format
        sortedCentroids - codebook of centroids to be used in future clustering
        subsetIdxes - return subset of indices; returns value only if
                      templateFlag is true
    """
    misc.vprint("Loading prior data.", verbose)
    misc.vprint("Please wait...", verbose)

    if not os.path.exists(priorVTKPath):
        raise IOError("Error: Prior data %s does not exist" % priorVTKPath)

    # Prior information

    priorTree = fibers.FiberTree()
    priorVTK, priorTree.no_of_fibers, priorTree.pts_per_fiber = \
        getFiberInfo(priorVTKPath)
    priorTree.convertFromVTK(priorVTK, priorTree.pts_per_fiber, verbose)

    # Get cluster labels + set number of fibers
    clusterCentroids, clusterArray = _getClusterInfo(priorVTK)

    # Get spatial information
    if templateFlag is True:
        subsetIdxes = _getSubset(clusterArray)

        centroidTree = priorTree.getFibers(subsetIdxes)
        centroidTree = fibers.convertFromTuple(centroidTree)
        _getScalarInfo(priorVTK, centroidTree, subsetIdxes,
                       centroidTree.pts_per_fiber, verbose)
        clusterArray = _addCentroidInfo(centroidTree, subsetIdxes,
                        clusterArray)

    else:
        centroidTree = priorTree.getFibers(range(priorTree.no_of_fibers))
        centroidTree = fibers.convertFromTuple(centroidTree)
        _getScalarInfo(priorVTK, centroidTree, range(priorTree.no_of_fibers),
                       centroidTree.pts_per_fiber, verbose)
        clusterArray = _addCentroidInfo(centroidTree,
                            range(priorTree.no_of_fibers), clusterArray)

        subsetIdxes = None

    misc.vprint("Finishined loading prior data.", verbose)

    del priorVTK, priorTree

    return centroidTree, clusterCentroids, clusterArray, subsetIdxes

def getFiberInfo(priorVTKPath):
    """
    Function to retrieve number of points from vtk polydata

    INPUT:
        priorVTKPath - path of .vtk polydata filer

    OUTPUT:
        no_of_fibers - number of fibers
        pts_per_fiber - number of samples along fiber
    """

    priorVTK = tractio.readVTK(priorVTKPath)
    no_of_fibers = priorVTK.GetNumberOfLines()
    pts_per_fiber = int(priorVTK.GetNumberOfPoints() / no_of_fibers)

    return priorVTK, no_of_fibers, pts_per_fiber

def _getSubset(clusterArray):
    """ *INTERNAL FUNCTION*
    Function to extract subset of fibers from each cluster. Used to
    subset template

    INPUT:
        clusterArray - array of cluster labels for each fiber

    OUTPUT:
        subsetIdxes - array of indices with subset from each fiber
    """

    subsetIdxes = []

    for cluster in np.unique(clusterArray):
        idx = np.where(clusterArray == cluster)[0]
        if len(idx) > 25:
            subsetIdx = np.random.choice(idx, 25, replace=False)
        else:
            subsetIdx = np.array(idx)

        subsetIdxes.extend(subsetIdx)

    return subsetIdxes

def _addCentroidInfo(centroidTree, subsetIdxes, clusterArray):
    """ *INTERNAL FUNCTION*
    Function to add centroid info to subset tree.

    INPUT:
        centroidTree - fiber tree containing subset fibers
        subsetIdxes - indices corresponding to subset of fibers
        clusterArray - array containing centroid info for all fibers

    OUTPUT:
        nClusterArray - array with new cluster info for subset
    """
    nClusterArray = clusterArray[subsetIdxes]

    for fidx in range(centroidTree.no_of_fibers):
            centroidTree.fiberTree[fidx][str(nClusterArray[fidx])] = \
                nClusterArray[fidx]

    return nClusterArray

def _getClusterInfo(priorVTK):
    """ *INTERNAL FUNCTION*
    Function to add cluster info to tree containing prior information. Also
    returns unique cluster centroids

    INPUT:
        priorVTK  - prior of .vtk polydata file

    OUTPUT:
        sortedCentroid - array of sorted centroids to be used for future clustering
    """

    # Get cluster information in array
    clusterLabels = priorVTK.GetCellData().GetArray('ClusterLabel')
    centroidLabels = priorVTK.GetCellData().GetArray('Centroid')
    clusterArray = numpy_support.vtk_to_numpy(clusterLabels)
    centroidArray = numpy_support.vtk_to_numpy(centroidLabels)

    for cluster in np.unique(clusterArray):
        idx = np.where(clusterArray == cluster)[0][0]

        # Sort centroids
        if cluster == 0:
            sortedCentroid = centroidArray[idx]
        else:
            sortedCentroid = np.vstack((sortedCentroid, centroidArray[idx]))

    del clusterLabels, centroidLabels, centroidArray

    return sortedCentroid, clusterArray

def _getScalarInfo(priorVTK, centroidTree, subsetIdx, pts_per_fiber=20,
                   verbose=0):
    """ *INTERNAL FUNCTION*
    Reads and converts scalar information stored in VTK to fiberTree

    INPUT:
        priorVTK - prior .vtk polydata file
        centroidTree - tree containing subset fiber data
        subsetIdx - subset of fibers to extract scalars from
        pts_per_fiber - number of points to sample along

    OUTPUT:
        none
    """

    for i in range(priorVTK.GetPointData().GetNumberOfArrays()):
        scalarType = priorVTK.GetPointData().GetArray(i).GetName()

        misc.vprint("Adding %s to fiber data" % scalarType, verbose)

        fidx = 0
        for idx in subsetIdx:
            if idx == 0:
                j = 0
            else:
                j = idx * pts_per_fiber

            for pidx in range(0, centroidTree.pts_per_fiber):
                centroidTree.fiberTree[fidx][pidx][scalarType] = \
                    priorVTK.GetPointData().GetArray(i).GetValue(j)
                j += 1

            fidx += 1

def loadEig(dirpath, eigvalFile, eigvecFile):
    """ WARNING: TO BE DEPRECATED
    Loads eigenvalue and eigenvector arrays from previously clustered data

    INPUT:
        dirpath - Directory path where eigenvalues & eigenvectors are stored
        eigvalArray - Array of eigenvalues to be loaded
        eigvecArray - Matrix of eigenvectos to be loaded

    OUTPUT:
        priorEigval - Numpy array of eigenvalues
        priorEigvec - Numpy array of eigenvectors
    """

    priorEigval = np.load(dirpath + '/' + eigvalFile)
    priorEigvec = np.load(dirpath + '/' + eigvecFile)

    return priorEigval, priorEigvec
