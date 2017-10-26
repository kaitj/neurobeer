""" prior.py

Module containing classes and functions pertaining to use of previously clustered data.

"""

import os
import numpy as np
import fibers, tractio
from vtk.util import numpy_support

def load(priorVTKPath, verbose=0):
    """
    Class used to load .vtk prior file.

    INPUT:
        priorVTKPath - absolute path to VTK file containing prior information to be used
        verbose - verbosity of function; defaults 0

    OUTPUT:
        priorTree - returns prior information stored in a tree format
        sortedCentroids - codebook of centroids to be used in future clustering
    """
    if verbose == 1:
        print '\nLoading prior data.'
        print '\nPlease wait...'

    if not os.path.exists(priorVTKPath):
        print 'Error: Prior data', priorVTKPath, 'does not exist.'
        raise IOError

    priorVTK = tractio.readVTK(priorVTKPath, verbose)

    priorPts = priorVTK.GetNumberOfPoints() / priorVTK.GetNumberOfLines()

    # Tree to store prior data
    priorTree = fibers.FiberTree()
    priorTree.convertFromVTK(priorVTK, pts_per_fiber=priorPts, verbose=0)

    # Get cluster labels
    priorCentroids = _getClusterInfo(priorVTK, priorTree)

    # Get scalar data
    _scalarFromVTK(priorVTK, priorTree, pts_per_fiber=priorPts, verbose=0)

    if verbose == 1:
        print '\nFinished loading prior data.'

    return priorTree, priorCentroids

def _getClusterInfo(priorVTK, priorTree):
    """ *INTERNAL FUNCTION*
    Function to add cluster info to tree containing prior information. Also returns unique
    cluster centroids

    INPUT:
        priorVTK  - prior of .vtk polydata file
        priorTree - tree containing prior fiber information

    OUTPUT:
        sortedCentroid - array of sorted centroids to be used for future clustering
    """

    # Get cluster information in array
    clusterLabels = priorVTK.GetCellData().GetArray('ClusterLabel')
    centroidLabels = priorVTK.GetCellData().GetArray('Centroid')
    clusterArray = numpy_support.vtk_to_numpy(clusterLabels)
    centroidArray = numpy_support.vtk_to_numpy(centroidLabels)
    del clusterLabels, centroidLabels

    # Add inforrmation to tree
    for fidx in range(priorTree.no_of_fibers):
        priorTree.fiberTree[fidx][str(clusterArray[fidx])] = clusterArray[fidx]

    # Find unique centroids
    for label in np.unique(clusterArray):
        if label == 0:
            idx = np.where(clusterArray == label)[0][0]
            priorTree.fiberTree['centroid'][label] = centroidArray[idx]
            sortedCentroid = centroidArray[idx]
        else:
            idx = np.where(clusterArray == label)[0][0]
            priorTree.fiberTree['centroid'][label] = centroidArray[idx]
            sortedCentroid = np.vstack((sortedCentroid, centroidArray[idx]))

    del clusterArray, centroidArray

    return sortedCentroid

def _scalarFromVTK(inputVTK, fiberTree, pts_per_fiber=20, verbose=0):
    """ *INTERNAL FUNCTION*
    Reads and converts scalar information stored in VTK to fiberTree

    INPUT:
        inputVTK - prior .vtk polydata file
        fiberTree - tree containing prior data
        pts_per_fiber - number of points to sample along fiber; defaults to 20
        verbose - verbosity of function; defaults 0

    OUTPUT:
        none
    """

    for i in range(inputVTK.GetPointData().GetNumberOfArrays()):
        scalarType = inputVTK.GetPointData().GetArray(i).GetName()
        idx = 0

        if verbose == 1:
            print "\nAdding %s to fiber data" % scalarType

        for fidx in range(0, fiberTree.no_of_fibers):
            for pidx in range(0, fiberTree.pts_per_fiber):
                fiberTree.fiberTree[fidx][pidx][scalarType] = \
                    inputVTK.GetPointData().GetArray(i).GetValue(idx)
                idx += 1

def loadEig(dirpath, eigvalFile, eigvecFile):
    """ *INTERNAL FUNCTION*
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
