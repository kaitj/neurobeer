""" distance.py

Module containing classes and functions related to calculating distance and
similarity measurements.

"""

import numpy as np
import scipy as sp

def _fiberDistance_internal(fiberMatrix, flip=False):
    """ *INTERNAL FUNCTION*
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers using MeanSquared method.

    INPUT:
        fiberMatrix1 - 3D matrix containing fiber spatial infomration
        fiberMatrix2 - 3D matrix containing fiber spatial information to be
                       compared with

    OUTPUT:
        distance - computed distance from single fiber to those within fiber
                    group
    """

    # Calculates the squared distance between fibers
    if flip is False:
        dx_sq = sp.spatial.distance.cdist(fiberMatrix[0, :, :],
            fiberMatrix[0, :, :], metric='sqeuclidean')
        dy_sq = sp.spatial.distance.cdist(fiberMatrix[1, :, :],
            fiberMatrix[1, :, :], metric='sqeuclidean')
        dz_sq = sp.spatial.distance.cdist(fiberMatrix[2, :, :],
            fiberMatrix[2, :, :], metric='sqeuclidean')
    else:
        dx_sq = sp.spatial.distance.cdist(np.fliplr(fiberMatrix[0, :, :]),
            fiberMatrix[0, :, :], metric='sqeuclidean')
        dy_sq = sp.spatial.distance.cdist(np.fliplr(fiberMatrix[1, :, :]),
            fiberMatrix[1, :, :], metric='sqeuclidean')
        dz_sq = sp.spatial.distance.cdist(np.fliplr(fiberMatrix[2, :, :]),
            fiberMatrix[2, :, :], metric='sqeuclidean')

    # Computed distance
    distance = np.sqrt(dx_sq + dy_sq + dz_sq)

    return distance

def _scalarDistance_internal(fiberScalarMatrix, flip=False):
    """ *INTERNAL FUNCTION*
    Computes the "distance" between the scalar values between one fiber and
    the fibers within a group (array) of fibers using MeanSquared method.

    INPUT:
        fiberScalarMatrix - array of scalar information pertaining to a group
                            of fibers
    OUTPUT:
        qDistance - computed scalar "distance" between fibers
    """

    # Calculates squared distance of scalars
    if flip is False:
        dq_sq = sp.spatial.distance.cdist(fiberScalarMatrix, fiberScalarMatrix,
            metric='sqeuclidean')
    else:
        dq_sq = sp.spatial.distance.cdist(np.fliplr(fiberScalarMatrix),
            fiberScalarMatrix, metric='sqeuclidean')

    qDistance = np.sqrt(dq_sq)

    return qDistance

def fiberDistance(fiberArray):
    """
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers. This function also handles equivalent fiber
    representations.

    INPUT:
        fiberMatrix - group of fibers for comparison

    OUTPUT:
        distance - minimum distance between group of fiber and single fiber
                            traversed in both directions
    """

    fiberMatrix = np.asarray(fiberArray)

    # Compute distances for fiber and fiber equivalent to fiber group
    distance1 = _fiberDistance_internal(fiberMatrix)
    distance2 = _fiberDistance_internal(fiberMatrix, flip=True)

    # Minimum distance more likely to be part of cluster; return distance
    distance = np.minimum(distance1, distance2)

    return distance

def scalarDistance(fiberScalarArray):
    """
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers. This function also handles equivalent fiber
    representations.

    INPUT:
        fiberScalarArray - array of scalar information pertaining to a group of fibers

    OUTPUT:
        distance - distance between group of fiber and single fiber

    TODO: Add functionality to calculate reverse fiber if necessary?
    """

    fiberScalarMatrix = np.array(fiberScalarArray)

    # Compute distances for fiber and fiber equivalent to fiber group
    distance1 = _scalarDistance_internal(fiberScalarMatrix)
    distance2 = _scalarDistance_internal(fiberScalarMatrix, flip=True)

    # Minimum distance more likely to be similar; return distance
    distance = np.minimum(distance1, distance2)

    return distance

def gausKernel_similarity(distance, sigmasq):
    """
    Computes the similarity using a Gaussian (RBF) kernel.

    INPUT:
        distance - Euclidean distance between points
        sigma - width of the kernel; adjust to alter sensitivity

    OUTPUT:
        similiarities - scalar values pertaining to the similarity of fiber
                                group of fibers (0 is dissimilar, 1 is identical)
    """

    # Computes similarity using a Gaussian kernel
    similarities = np.exp(-distance / sigmasq)

    return similarities
