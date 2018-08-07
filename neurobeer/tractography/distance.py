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
        flip - flag to flip fiber

    OUTPUT:
        distance - Matrix containing distance between fibers
    """

    # Comparison of single fiber group
    if fiberMatrix2 is None:
        distance = np.empty((fiberMatrix1.shape[1], fiberMatrix1.shape[1]),
                             dtype=np.float32)

        # Calculates the squared distance between fibers
        if flip is False:
            for i in range(fiberMatrix1.shape[1]):
                distance[i] = np.sqrt((np.linalg.norm(fiberMatrix1[0, i, :] -
                                                      fiberMatrix1[0, :, :],
                                                      axis=1))**2 +
                                      (np.linalg.norm(fiberMatrix1[1, i, :] -
                                                      fiberMatrix1[1, :, :],
                                                      axis=1))**2 +
                                      (np.linalg.norm(fiberMatrix1[2, i, :] -
                                                      fiberMatrix1[2, :, :],
                                                      axis=1))**2)
        # Flipped fiber
        else:
            for i in range(fiberMatrix1.shape[1]):
                distance[i] = np.sqrt((np.linalg.norm(
                                       np.flip(fiberMatrix1[0, i, :], axis=0) -
                                               fiberMatrix1[0, :, :],
                                               axis=1))**2 +
                                      (np.linalg.norm(
                                       np.flip(fiberMatrix1[1, i, :], axis=0) -
                                               fiberMatrix1[1, :, :],
                                               axis=1))**2 +
                                      (np.linalg.norm(
                                       np.flip(fiberMatrix1[2, i, :], axis=0) -
                                               fiberMatrix1[2, :, :],
                                               axis=1))**2)
        del fiberMatrix1
    # Comparison between two fiber groups
    else:
        distance = np.empty((fiberMatrix1.shape[1], fiberMatrix2.shape[1]),
                             dtype=np.float32)

        # Calculates the squared distance between fibers
        if flip is False:
            for i in range(fiberMatrix1.shape[1]):
                distance[i] = np.sqrt((np.linalg.norm(fiberMatrix1[0, i, :] -
                                                      fiberMatrix2[0, :, :],
                                                      axis=1))**2 +
                                      (np.linalg.norm(fiberMatrix1[1, i, :] -
                                                      fiberMatrix2[1, :, :],
                                                      axis=1))**2 +
                                      (np.linalg.norm(fiberMatrix1[2, i, :] -
                                                      fiberMatrix2[2, :, :],
                                                      axis=1))**2)
        # Flipped fiber
        else:
            for i in range(fiberMatrix1.shape[1]):
                distance[i] = np.sqrt((np.linalg.norm(
                                       np.flip(fiberMatrix1[0, i, :], axis=0) -
                                               fiberMatrix2[0, :, :],
                                               axis=1))**2 +
                                      (np.linalg.norm(
                                       np.flip(fiberMatrix1[1, i, :], axis=0) -
                                               fiberMatrix2[1, :, :],
                                               axis=1))**2 +
                                      (np.linalg.norm(
                                       np.flip(fiberMatrix1[2, i, :], axis=0) -
                                               fiberMatrix2[2, :, :],
                                               axis=1))**2)
        del fiberMatrix1, fiberMatrix2

    return distance

def _scalarDistance_internal(fiberScalarMatrix1, fiberScalarMatrix2=None,
                             flip=False):
    """ *INTERNAL FUNCTION*
    Computes the "distance" between the scalar values between one fiber and
    the fibers within a group (array) of fibers.

    INPUT:
        fiberScalarMatrix1 - array of scalar information pertaining to a group
                             of fibers
        fiberScalarMatrix2 - array of scalar information pertaining to a group
                             of fibers for comparison
        flip - flag to flip fiber
    OUTPUT:
        qDistance - computed scalar "distance" between fibers
    """

    # Calculates squared distance of scalars
    if fiberScalarMatrix2 is None:
        qDistance = np.empty((fiberScalarMatrix1.shape[1],
                              fiberScalarMatrix1.shape[1]), dtype=np.float32)

        # Calculates the squared distance between fibers
        if flip is False:
            for i in range(fiberScalarMatrix1.shape[1]):
                qDistance[i] = np.sqrt((np.linalg.norm(
                                            fiberScalarMatrix1[0, i, :] -
                                            fiberScalarMatrix1[0, :, :],
                                            axis=1))**2 +
                                       (np.linalg.norm(
                                            fiberScalarMatrix1[1, i, :] -
                                            fiberScalarMatrix1[1, :, :],
                                            axis=1))**2 +
                                       (np.linalg.norm(
                                            fiberScalarMatrix1[2, i, :] -
                                            fiberScalarMatrix1[2, :, :],
                                            axis=1))**2)
        # Flip fiber
        else:
            for i in range(fiberScalarMatrix1.shape[1]):
                qDistance[i] = np.sqrt((np.linalg.norm(
                                        np.flip(fiberScalarMatrix1[0, i, :],
                                                axis=0) -
                                                fiberScalarMatrix1[0, :, :],
                                                axis=1))**2 +
                                       (np.linalg.norm(
                                        np.flip(fiberScalarMatrix1[1, i, :],
                                                axis=0) -
                                                fiberScalarMatrix1[1, :, :],
                                                axis=1))**2 +
                                       (np.linalg.norm(
                                        np.flip(fiberScalarMatrix1[2, i, :],
                                                axis=0) -
                                                fiberScalarMatrix1[2, :, :],
                                                axis=1))**2)
        del fiberScalarMatrix1
    # Comparison between two fiber groups
    else:
        distance = np.empty((fiberScalarMatrix1.shape[1],
                             fiberScalarMatrix2.shape[1]), dtype=np.float32)

        # Calculates the squared distance between fibers
        if flip is False:
            for i in range(fiberScalarMatrix1.shape[1]):
                qDistance[i] = np.sqrt((np.linalg.norm(
                                            fiberScalarMatrix1[0, i, :] -
                                            fiberScalarMatrix2[0, :, :],
                                            axis=1))**2 +
                                       (np.linalg.norm(
                                            fiberScalarMatrix1[1, i, :] -
                                            fiberScalarMatrix2[1, :, :],
                                            axis=1))**2 +
                                       (np.linalg.norm(
                                            fiberScalarMatrix1[2, i, :] -
                                            fiberScalarMatrix2[2, :, :],
                                            axis=1))**2)
        # Flip fiber
        else:
            for i in range(fiberScalarMatrix1.shape[1]):
                qDistance[i] = np.sqrt((np.linalg.norm(
                                        np.flip(fiberScalarMatrix1[0, i, :],
                                                axis=0) -
                                                fiberScalarMatrix2[0, :, :],
                                                axis=1))**2 +
                                       (np.linalg.norm(
                                        np.flip(fiberScalarMatrix1[1, i, :],
                                                axis=0) -
                                                fiberScalarMatrix2[1, :, :],
                                                axis=1))**2 +
                                       (np.linalg.norm(
                                        np.flip(fiberScalarMatrix1[2, i, :],
                                                axis=0) -
                                                fiberScalarMatrix2[2, :, :],
                                                axis=1))**2)

        del fiberScalarMatrix1, fiberScalarMatrix2

    return distance

def fiberDistance(fiberArray1, fiberArray2=None):
    """
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers. This function also handles equivalent fiber
    representations.

    INPUT:
        fiberArray1 - group of fibers for comparison
        fiberArray2 - group of fibers to compare fiberArray1 to, if applicable

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

def scalarDistance(fiberScalarArray1, fiberScalarArray2=None):
    """
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers. This function also handles equivalent fiber
    representations.

    INPUT:
        fiberScalarArray1 - array of scalar information pertaining to a group of
                            fibers
        fiberScalarArray2 - array of scalar information pertaining to a group of
                            fibers if applicable

    OUTPUT:
        distance - distance between group of fiber and single fiber
    """
    if fiberScalarArray2 is None:
        fiberScalarArray1 = np.asarray(fiberScalarArray1, dtype=np.float32)

        # Compute distances for fiber and fiber equivalent to fiber group
        distance1 = _scalarDistance_internal(fiberScalarArray1)
        distance2 = _scalarDistance_internal(fiberScalarArray1, flip=True)
        del fiberScalarArray1

    else:
        fiberScalarArray1 = np.asarray(fiberScalarArray1, dtype=np.float32)
        fiberScalarArray2 = np.asarray(fiberScalarArray2, dtype=np.float32)

        # Compute distances between two fiber groups
        distance1 = _scalarDistance_internal(fiberScalarArray1,
                                             fiberScalarArray2)
        distance2 = _scalarDistance_internal(fiberScalarArray1,
                                             fiberScalarArray2,
                                             flip=True)
        del fiberScalarArray1, fiberScalarArray2

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
