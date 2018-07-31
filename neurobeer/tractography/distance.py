""" distance.py

Module containing classes and functions related to calculating distance and
similarity measurements.

"""

import numpy as np
import scipy as sp

def _fiberDistance_internal(fiberMatrix1, fiberMatrix2=None, flip=False):
    """ *INTERNAL FUNCTION*
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers using MeanSquared method.

    INPUT:
        fiberMatrix1 - 3D matrix containing fiber spatial infomration
        fiberMatrix2 - 3D matrix containing fiber spatial information for
                       comparison

    OUTPUT:
        distance - Matrix containing distance between fibers
    """

    # Initialize array
    if fiberMatrix2 is None:
        distance = np.zeros((fiberMatrix1.shape[1], fiberMatrix1.shape[1]),
                    dtype=np.float32)

        # Calculates the squared distance between fibers
        if flip is False:
            for i in range(3):
                distance += sp.spatial.distance.cdist(fiberMatrix1[i, :, :],
                    fiberMatrix1[i, :, :], metric='sqeuclidean')
        else:
            for i in range(3):
                distance += sp.spatial.distance.cdist(np.fliplr(fiberMatrix1[i, :, :]),
                    fiberMatrix1[i, :, :], metric='sqeuclidean')
    else:
        distance = np.zeros((fiberMatrix1.shape[1], fiberMatrix2.shape[1]),
                    dtype=np.float32)

        # Calculates the squared distance between fibers
        if flip is False:
            for i in range(3):
                distance += sp.spatial.distance.cdist(fiberMatrix1[i, :, :],
                    fiberMatrix2[i, :, :], metric='sqeuclidean')
        else:
            for i in range(3):
                distance += sp.spatial.distance.cdist(np.fliplr(fiberMatrix1[i, :, :]),
                    fiberMatrix2[i, :, :], metric='sqeuclidean')

    # Computed distance
    distance = np.sqrt(distance)

    # Delete variables no longer needed
    del fiberMatrix1, fiberMatrix2

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
    qDistance = np.zeros((fiberScalarMatrix.shape[1],
                fiberScalarMatrix.shape[1]), dtype=np.float32)

    if flip is False:
        qDistance += sp.spatial.distance.cdist(fiberScalarMatrix, fiberScalarMatrix,
            metric='sqeuclidean')
    else:
        qDistance += sp.spatial.distance.cdist(np.fliplr(fiberScalarMatrix),
            fiberScalarMatrix, metric='sqeuclidean')

    qDistance = np.sqrt(qDistance)

    # Delete variables no longer needed
    del fiberScalarMatrix

    return qDistance

def fiberDistance(fiberArray1, fiberArray2=None):
    """ fiberArray
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers. This function also handles equivalent fiber
    representations.

    INPUT:
        fiberMatrix - group of fibers for comparison

    OUTPUT:
        distance - minimum distance between group of fiber and single fiber
                   traversed in both directions
    """

    if fiberArray2 is None:
        fiberArray1 = np.asarray(fiberArray1, dtype=np.float32)

        # Compute distances for fiber and flipped fiber of group
        distance1 = _fiberDistance_internal(fiberArray1)
        distance2 = _fiberDistance_internal(fiberArray1, flip=True)

        del fiberArray1
    else:
        fiberArray1 = np.asarray(fiberArray1, dtype=np.float32)
        fiberArray2 = np.asarray(fiberArray2, dtype=np.float32)

        # Compute distances between two fiber groups
        distance1 = _fiberDistance_internal(fiberArray1, fiberArray2)
        distance2 = _fiberDistance_internal(fiberArray1, fiberArray2, flip=True)

        del fiberArray1, fiberArray2

    # Minimum distance more likely to be part of cluster; return distance
    distance = np.minimum(distance1, distance2)

    del distance1, distance2

    return distance

def scalarDistance(fiberScalarArray):
    """
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers. This function also handles equivalent fiber
    representations.

    INPUT:
        fiberScalarArray - array of scalar information pertaining to a group of
                           fibers

    OUTPUT:
        distance - distance between group of fiber and single fiber

    TODO: Add functionality to calculate reverse fiber if necessary?
    """

    fiberScalarArray = np.array(fiberScalarArray, dtype=np.float32)

    # Compute distances for fiber and fiber equivalent to fiber group
    distance1 = _scalarDistance_internal(fiberScalarArray)
    distance2 = _scalarDistance_internal(fiberScalarArray, flip=True)

    # Minimum distance more likely to be similar; return distance
    distance = np.minimum(distance1, distance2)

    del distance1, distance2, fiberScalarArray

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

    del distance, sigmasq

    return similarities
