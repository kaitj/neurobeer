""" cluster.py

Module containing classes and functions used to cluster fibers and modify
parameters pertaining to clusters.

"""

import numpy as np
import fibers, distance, scalars
from joblib import Parallel, delayed

class Cluster:
    """ Clustering of whole-brain tractography data from subject """

    def cluster(inputVTK, k_clusters, sigma, no_of_jobs):
        """
        Clustering of fibers based on pairwise fiber similarity

            INPUT:
                inputVTK - input polydata file
                k_clusters - number of clusters via k-means clustering
                sigma - width of kernel; adjust to alter sensitivity
                no_of_jobs - processes to use to perform computation
        """

        noFibers = inputVTK.GetNumberOfLines()
        if noFibers == 0:
            print "ERROR: Input data has 0 fibers!"
            return
        else:
            print "Starting clustering..."
            print "No. of fibers:", noFibers
            print "No. of clusters:", k_clusters


def _pairwiseDistance_matrix(inputVTK, sigma, no_of_jobs):
    """ An internal function used to compute an NxN distance matrix for all
    fibers (N) in the input data

    INPUT:
        inputVTK - input polydata file
        sigma - width of kernel; adjust to alter sensitivity
        no_of_jobs - processes to use to perform computation

    OUTPUT:
        distances - NxN matrix containing distances between fibers
    """

    fiberArray = fibers.FiberArray()
    fiberArray.convertFromVTK(inputVTK, pts_per_fiber=20)

    distances = Parallel(n_jobs=no_of_jobs, verbose=0)(
            delayed(distance.fiberDistance)(fiberArray.getFiber(fidx),
                    fiberArray)
            for fidx in range(0, fiberArray.no_of_fibers))

    distances = np.array(distances)

    return distances

def _pairwiseSimilarity_matrix(inputVTK, sigma, no_of_jobs):
    """ An internal function used to compute an NxN similarity matrix for all
    fibers (N) in the input data.

    INPUT:
        inputVTK - input polydata file
        sigma - width of kernel; adjust to alter sensitivity
        no_of_jobs - proccesses to use to perform computation

    OUTPUT:
        similarity - NxN matrix containing similarity between fibers
    """

    distances = _pairwiseDistance_matrix(inputVTK, sigma, no_of_jobs)

    sigmasq = np.square(sigma)
    similarities = distance.gausKernel_similarity(distances, sigmasq)

    return similarities

def _pairwiseQSimilarity_matrix(inputVTK, scalarData, scalarType, no_of_jobs):
    """ An internal function used to compute the cross-correlation between
    quantitative metrics along a fiber.

    INPUT:
        inputVTK - input polydata file
        no_of_jobs - processes to use to perform computation

    OUTPUT:
        qSimilarity - NxN matrix containing correlation values between fibers
    """

    fiberArray = fibers.FiberArray()
    scalarArray = scalars.FiberArrayScalar()

    fiberArray.convertFromVTK(inputVTK, no_of_pts=20)
    scalarArray.addScalar(inputVTK, fiberArray, scalarData, scalarType)

    qSimilarity = Parallel(n_jobs=no_of_jobs, verbose=0)(
            delayed(np.corrcoef[0][1])(
                    scalarArray.getScalar(fiberArray, scalarType)[fidx],
                    scalarArray)
            for fidx in range(0, fiberArray.no_of_fibers))

    qSimilarity = np.array(qSimilarity)

    return qSimilarity
