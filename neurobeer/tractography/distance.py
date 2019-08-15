""" distance.py

Module containing classes and functions related to calculating distance and
similarity measurements.

"""

import numpy as np
from joblib import Parallel, delayed
from joblib.pool import has_shareable_memory

def _calcDistance(fiberMatrix1, fiberMatrix2):
    """ *INTERNAL FUNCTION*
    Computes average Euclidean distance

    INPUT:
        fiberMatrix1 - 3D matrix containing fiber spatial infomration
        fiberMatrix2 - 3D matrix containing fiber spatial information for
                       comparison

    OUTPUT:
        Average Euclidean distance of sample points
    """
    return np.mean(np.linalg.norm(np.subtract(fiberMatrix1, fiberMatrix2),
            axis=0), axis=1)

def _calcQDistance(fiberMatrix1, fiberMatrix2):
    """ *INTERNAL FUNCTION*
    Computes average Euclidean distance

    INPUT:
        fiberMatrix1 - 3D matrix containing fiber quantitative infomration
        fiberMatrix2 - 3D matrix containing fiber quantitative information for
                       comparison

    OUTPUT:
        Average "Euclidean" distance of quantitative values
    """
    return np.mean(np.linalg.norm(fiberMatrix1, fiberMatrix2), axis=1)

def _fiberDistance_internal(fiberMatrix1, fiberMatrix2, flip=False, n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers.

    INPUT:
        fiberMatrix1 - 3D matrix containing fiber spatial infomration
        fiberMatrix2 - 3D matrix containing fiber spatial information for
                       comparison
        flip - flag to flip fiber
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        distance - Matrix containing distance between fibers
    """

    # Calculates the avg Euclidean distance between fibers
    if flip is False:
        distance = Parallel(n_jobs=n_jobs, backend='threading')(
            delayed(_calcDistance, has_shareable_memory)(
                fiberMatrix1[:, i, None], fiberMatrix2)
            for i in range(fiberMatrix1.shape[1]))

    # Flipped fiber
    else:
        distance = Parallel(n_jobs=n_jobs, backend='threading')(
            delayed(_calcDistance, has_shareable_memory)(
                np.flip(fiberMatrix1[:, i, None], axis=2), fiberMatrix2)
            for i in range(fiberMatrix1.shape[1]))

    return distance

def _scalarDistance_internal(fiberScalarMatrix1, fiberScalarMatrix2,
                             flip=False, n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Computes the "distance" between the scalar values between one fiber and
    the fibers within a group (array) of fibers.

    INPUT:
        fiberScalarMatrix1 - array of scalar information pertaining to a group
                             of fibers
        fiberScalarMatrix2 - array of scalar information pertaining to a group
                             of fibers for comparison
        flip - flag to flip fiber
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        qDistance - computed scalar "distance" between fibers
    """

    # Calculates squared distance of scalars
    qDistance = np.empty((fiberScalarMatrix1.shape[0],
                         fiberScalarMatrix2.shape[0]), dtype=np.float32)

    # Calculates the mean distance between fiber metrics
    if flip is False:
        qDistance = Parallel(n_jobs=n_jobs, backend='threading')(
            delayed(_calcQDistance, has_shareable_memory)(
                fiberScalarMatrix1[i, :], fiberScalarMatrix2)
            for i in range(fiberScalarMatrix1.shape[0]))

    # Flip fiber
    else:
        qDistance = Parallel(n_jobs=n_jobs, backend='threading')(
            delayed(_calcQDistance, has_shareable_memory)(
                np.flip(fiberScalarMatrix1[i, :], fiberScalarMatrix2))
            for i in range(fiberScalarMatrix1.shape[0]))

    return qDistance

def fiberDistance(fiberArray1, fiberArray2=None, n_jobs=-1):
    """
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers. This function also handles equivalent fiber
    representations.

    INPUT:
        fiberArray1 - group of fibers for comparison
        fiberArray2 - group of fibers to compare fiberArray1 to, if applicable
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        distance - minimum distance between group of fiber and single fiber
                   traversed in both directions
    """

    if fiberArray2 is None:
        fiberArray1 = np.asarray(fiberArray1, dtype=np.float32)

        # Compute distances for fiber and flipped fiber of group
        distance1 = _fiberDistance_internal(fiberArray1, fiberArray1,
                                            n_jobs=n_jobs)
        distance2 = _fiberDistance_internal(fiberArray1, fiberArray1, flip=True,
                                            n_jobs=n_jobs)
        del fiberArray1

    else:
        fiberArray1 = np.asarray(fiberArray1, dtype=np.float32)
        fiberArray2 = np.asarray(fiberArray2, dtype=np.float32)

        # Compute distances between two fiber groups
        distance1 = _fiberDistance_internal(fiberArray1, fiberArray2,
                                            n_jobs=n_jobs)
        distance2 = _fiberDistance_internal(fiberArray1, fiberArray2, flip=True,
                                            n_jobs=n_jobs)
        del fiberArray1, fiberArray2

    # Minimum distance more likely to be part of cluster; return distance
    distance = np.minimum(distance1, distance2)
    del distance1, distance2

    return distance

def scalarDistance(fiberScalarArray1, fiberScalarArray2=None, n_jobs=-1):
    """
    Computes the distance between one fiber and individual fibers within a
    group (array) of fibers. This function also handles equivalent fiber
    representations.

    INPUT:
        fiberScalarArray1 - array of scalar information pertaining to a group of
                            fibers
        fiberScalarArray2 - array of scalar information pertaining to a group of
                            fibers if applicable
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        distance - distance between group of fiber and single fiber
    """
    if fiberScalarArray2 is None:
        fiberScalarArray1 = np.asarray(fiberScalarArray1, dtype=np.float32)

        # Compute distances for fiber and fiber equivalent to fiber group
        distance1 = _scalarDistance_internal(fiberScalarArray1,
                                             fiberScalarArray1, n_jobs=n_jobs)
        distance2 = _scalarDistance_internal(fiberScalarArray1,
                                             fiberScalarArray1,
                                             flip=True, n_jobs=n_jobs)
        del fiberScalarArray1

    else:
        fiberScalarArray1 = np.asarray(fiberScalarArray1, dtype=np.float32)
        fiberScalarArray2 = np.asarray(fiberScalarArray2, dtype=np.float32)

        # Compute distances between two fiber groups
        distance1 = _scalarDistance_internal(fiberScalarArray1,
                                             fiberScalarArray2, n_jobs=n_jobs)
        distance2 = _scalarDistance_internal(fiberScalarArray1,
                                             fiberScalarArray2,
                                             flip=True, n_jobs=n_jobs)
        del fiberScalarArray1, fiberScalarArray2

    # Minimum distance more likely to be similar; return distance
    distance = np.minimum(distance1, distance2)
    del distance1, distance2

    return distance

def gausKernel_similarity(distance, sigma):
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
    similarities = np.exp(-np.square(distance) / np.square(sigma))
    del distance, sigma

    return similarities
