""" distance.py

Module containing classes and functions related to calculating distance and
similarity measurements.

"""

import numpy as np

def _fiberDistance_internal(fiberMatrix1, fiberMatrix2=None, flip=False):
    """ *INTERNAL FUNCTION*
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers.

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
            for i in range(fiberMatrix1.shape[1]):
                for j in range(fiberMatrix1.shape[1]):
                    if j <= i:
                        continue
                    distance[i, j] = np.linalg.norm(fiberMatrix1[:, i, :] -
                                                    fiberMatrix1[:, j, :])

            ind_lower = np.tril_indices(fiberMatrix1.shape[1], -1)
            distance[ind_lower] = distance.T[ind_lower]

        else:
            for i in range(fiberMatrix1.shape[1]):
                for j in range(fiberMatrix1.shape[1]):
                    if j <= i:
                        continue
                    distance[i, j] = np.linalg.norm(np.fliplr(
                                                    fiberMatrix1[:, i, :]) -
                                                    fiberMatrix1[:, j, :])

            ind_lower = np.tril_indices(fiberMatrix1.shape[1], -1)
            distance[ind_lower] = distance.T[ind_lower]

        del fiberMatrix1
    # Comparison between two fiber groups
    else:
        distance = np.zeros((fiberMatrix1.shape[1], fiberMatrix2.shape[1]),
                             dtype=np.float32)

        # Calculates the squared distance between fibers
        if flip is False:
            for i in range(fiberMatrix1.shape[1]):
                for j in range(fiberMatrix2.shape[1]):
                    distance[i, j] = np.linalg.norm(fiberMatrix1[:, i, :] -
                                                    fiberMatrix2[:, j, :])

        else:
            for i in range(fiberMatrix1.shape[1]):
                for j in range(fiberMatrix2.shape[1]):
                    distance[i, j] = np.linalg.norm(np.fliplr(
                                                fiberMatrix1[:, i, :]) -
                                                fiberMatrix1[:, j, :])

        del fiberMatrix1, fiberMatrix2

    return distance

def _scalarDistance_internal(fiberScalarMatrix1, fiberScalarMatrix2=None,
                             flip=False):
    """ *INTERNAL FUNCTION*
    Computes the "distance" between the scalar values between one fiber and
    the fibers within a group (array) of fibers.

    INPUT:
        fiberScalarMatrix - array of scalar information pertaining to a group
                            of fibers
    OUTPUT:
        qDistance - computed scalar "distance" between fibers
    """

    # Calculates squared distance of scalars
    if fiberScalarMatrix2 is None:
        qDistance = np.zeros((fiberScalarMatrix1.shape[1],
                             fiberScalarMatrix2.shape[1]), dtype=np.float32)

        # Calculates the squared distance between fibers
        if flip is False:
            for i in range(fiberScalarMatrix1.shape[1]):
                for j in range(fiberScalarMatrix1.shape[1]):
                    if j <= i:
                        continue
                    qDistance[i, j] = np.linalg.norm(
                                            fiberScalarMatrix1[:, i, :] -
                                            fiberScalarMatrix1[:, j, :])

            ind_lower = np.tril_indices(fiberScalarMatrix1.shape[1], -1)
            qDistance[ind_lower] = qDistance.T[ind_lower]

        else:
            for i in range(fiberScalarMatrix1.shape[1]):
                for j in range(fiberScalarMatrix2.shape[1]):
                    if j <= i:
                        continue
                    qDistance[i, j] = np.linalg.norm(np.fliplr(
                                                fiberScalarMatrix1[:, i, :]) -
                                                fiberScalarMatrix1[:, j, :])

            ind_lower = np.tril_indices(fiberScalarMatrix1.shape[1], -1)
            qDistance[ind_lower] = qDistance.T[ind_lower]

        del fiberScalarMatrix1
    # Comparison between two fiber groups
    else:
        distance = np.zeros((fiberScalarMatrix1.shape[1],
                             fiberScalarMatrix2.shape[1]), dtype=np.float32)

        # Calculates the squared distance between fibers
        if flip is False:
            for i in range(fiberScalarMatrix1.shape[1]):
                for j in range(fiberScalarMatrix2.shape[1]):
                    distance[i, j] = np.linalg.norm(
                                        fiberScalarMatrix1[:, i, :] -
                                        fiberScalarMatrix2[:, j, :])

        else:
            for i in range(fiberScalarMatrix1.shape[1]):
                for j in range(fiberScalarMatrix2.shape[1]):
                    distance[i, j] = np.linalg.norm(np.fliplr(
                                                fiberScalarMatrix1[:, i, :]) -
                                                fiberScalarMatrix2[:, j, :])

        del fiberScalarMatrix1, fiberScalarMatrix2

    return distance

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
        fiberArray1 = np.asarray(fiberArray1)

        # Compute distances for fiber and flipped fiber of group
        distance1 = _fiberDistance_internal(fiberArray1)
        distance2 = _fiberDistance_internal(fiberArray1, flip=True)
        del fiberArray1

    else:
        fiberArray1 = np.asarray(fiberArray1)
        fiberArray2 = np.asarray(fiberArray2)

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

    fiberScalarArray = np.array(fiberScalarArray)

    # Compute distances for fiber and fiber equivalent to fiber group
    distance1 = _scalarDistance_internal(fiberScalarArray)
    distance2 = _scalarDistance_internal(fiberScalarArray, flip=True)
    del fiberScalarArray

    # Minimum distance more likely to be similar; return distance
    distance = np.minimum(distance1, distance2)
    del distance1, distance2

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
