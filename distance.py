""" distance.py

Module containing classes and functions related to calculating distance and
similarity measurements.

"""

import numpy as np

def _fiberDistance_internal(fiber, fiberArray):
    """ Computes the distance between one fiber and individual fibers within a
    group (array) of fibers using MeanSquared method.

    INPUT:
        fiber - single fiber to be compared
        fiberArray - group of fibers that lone fiber is to be compared to

    OUTPUT:
        distance - computed distance from single fiber to those within fiber
        group
    """

    # Calculates the distance between points of each fiber
    dx = fiberArray.fiberArray_x - fiber.x
    dy = fiberArray.fiberArray_y - fiber.y
    dz = fiberArray.fiberArray_z - fiber.z

    # Squares calculated distances
    dx_sq = np.square(dx)
    dy_sq = np.square(dy)
    dz_sq = np.square(dz)

    # Sums the squared values of the distances in each axis
    distance = dx_sq + dy_sq + dz_sq

    # Sum along fiber and normalize
    distance = np.sum(distance, 1)
    distance = distance / float(fiberArray.pts_per_fiber)

    return distance

def fiberDistance(fiber, fiberArray):
    """ Computes the distance between one fiber and individual fibers within a
    group (array) of fibers. This function also handles equivalent fiber
    representations. 

    INPUT:
        fiber - single fiber to be compared
        fiberArray - group of fibers that lone fiber is to be compared to

    OUTPUT:
        distance - minimum distance between group of fiber and single fiber
        traversed in both directions
    """

    # Get same fiber in reverse order
    fiberEquiv = fiber.getReverseFiber()

    # Compute distances for fiber and fiber equivalent to fiber group
    distance1 = _fiberDistance_internal(fiber, fiberArray)
    distance2 = _fiberDistance_internal(fiberEquiv, fiberArray)

    # Minimum distance more likely to be part of cluster; return distance
    distance = np.minimum(distance1, distance2)

    return distance

def gausKernel_similarity(distance, sigmasq=100):
    """ Computes the similarity using a Gaussian kernel

    INPUT:
        distance - Euclidean distance between points
        sigma - width of the kernel; adjust to alter sensitivity

    OUTPUT:
        similiarities - scalar values pertaining to the similarity of fiber to
        group of fibers (0 is dissimilar, 1 is identical)
    """

    # Computes similarity using a Gaussian kernel
    similarities = np.exp(-distance / (sigmasq))

    return similarities
