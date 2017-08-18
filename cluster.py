""" cluster.py

Module containing classes and functions used to cluster fibers and modify
parameters pertaining to clusters.

"""

import numpy as np
import scipy.cluster
from joblib import Parallel, delayed
import fibers, distance, scalars
import vtk

import sklearn.cluster, sklearn.preprocessing

def spectralClustering(inputVTK, scalarDataList=[], scalarTypeList=[], scalarWeightList=[],
                                    k_clusters=4, no_of_eigvec=20, sigma=0.4, no_of_jobs=2):
        """
        Clustering of fibers based on pairwise fiber similarity
        See paper: "A tutorial on spectral clustering" (von Luxburg, 2007)

            INPUT:
                inputVTK - input polydata file
                k_clusters - number of clusters via k-means clustering
                sigma - width of kernel; adjust to alter sensitivity
                no_of_jobs - processes to use to perform computation

        TODO: weighted quantitative clustering
        """

        if no_of_eigvec == 1:
            print "\nClustering cannot be performed with single eigenvector!"
            return

        noFibers = inputVTK.GetNumberOfLines()
        if noFibers == 0:
            print "\nERROR: Input data has 0 fibers!"
            return
        else:
            print "\nStarting clustering..."
            print "No. of fibers:", noFibers
            print "No. of clusters:", k_clusters

        # 1. Compute similarty matrix
        W = _weightedSimilarity(inputVTK, scalarDataList, scalarTypeList, scalarWeightList,
                                                    sigma, no_of_jobs)

        # 2. Compute degree matrix
        D = _degreeMatrix(W)

        # 3. Compute unnormalized Laplacian
        L = D - W

        # 4. Compute normalized Laplacian (random-walk)
        Lrw = np.dot(np.diag(np.divide(1, np.sum(D, 0))), L)

        # 5. Compute eigenvalues and eigenvectors of generalized eigenproblem
        # vectors are columns ([:, n]) of matrix
        eigval, eigvec = np.linalg.eig(Lrw)

        # 6. Compute information for clustering using "N" number of smallest eigenvalues
        # Skip first eigenvector, no information provided for clustering???
        U = eigvec[:, 0:no_of_eigvec]
        U = U.astype('float')

        # 7. Find clusters using K-means clustering
        # Sort centroids by eigenvector order
        centroids, clusterIdx = scipy.cluster.vq.kmeans2(U, k_clusters, minit='points')

        if no_of_eigvec == 1:
            print('Not enough eigenvectors selected!')
        elif no_of_eigvec == 2:
            colour = _cluster_to_rgb(U)
        else:
            colour = _cluster_to_rgb(centroids)

        # 8. Return results
        outputData = inputVTK
        outputPolydata = _format_outputVTK(outputData, clusterIdx, colour, U)

        # 9. Also add measurements from those used to cluster
        for i in range(len(scalarDataList)):
            outputPolydata = addScalarToVTK(outputPolydata, scalarDataList[i], scalarTypeList[i])

        return outputPolydata, clusterIdx, colour, centroids

def addScalarToVTK(polyData, scalarData, scalarType):
    """ Add scalar to polydata points to be converted to .vtk file.

    This function is different from scalars.addScalar, which only considers point
    used in sampling of fiber.
    """

    data = vtk.vtkFloatArray()
    data.SetName(scalarType.split('/', -1)[-1])

    for pidx in range(0, polyData.GetNumberofPoints()):
        data.InsertNextValue(float(scalarData[pidx]))

    polyData.GetPointData().AddArray(data)

    return polyData

def _pairwiseDistance_matrix(inputVTK, no_of_jobs):
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
    # Normalize between 0 and 1
    distances = sklearn.preprocessing.MinMaxScaler().fit_transform(distances)

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

    distances = _pairwiseDistance_matrix(inputVTK, no_of_jobs)

    sigmasq = np.square(sigma)
    similarities = distance.gausKernel_similarity(distances, sigmasq)

    similarities = np.array(similarities)

    return similarities

def _pairwiseQDistance_matrix(inputVTK, scalarData, scalarType, no_of_jobs):
    """ An internal function used to compute the "distance" between quantitative
    points along a fiber. """

    fiberArray = fibers.FiberArray()
    fiberArray.convertFromVTK(inputVTK, pts_per_fiber=20)
    no_of_fibers = fiberArray.no_of_fibers
    scalarArray = scalars.FiberArrayScalar()
    scalarArray.addScalar(inputVTK, fiberArray, scalarData, scalarType)

    qDistances = Parallel(n_jobs=no_of_jobs, verbose=0)(
            delayed(distance.scalarDistance)(
                scalarArray.getScalar(fiberArray, fidx, scalarType),
                scalarArray.getScalars(fiberArray, range(no_of_fibers), scalarType))
            for fidx in range(0, no_of_fibers)
    )

    qDistances = np.array(qDistances)

    # Normalize distance measurements
    qDistances = sklearn.preprocessing.MinMaxScaler().fit_transform(qDistances)

    return qDistances

def _pairwiseQSimilarity_matrix(inputVTK, scalarData, scalarType, sigma, no_of_jobs):
    """ An internal function used to compute the cross-correlation between
    quantitative metrics along a fiber.

    INPUT:
        inputVTK - input polydata file
        no_of_jobs - processes to use to perform computation

    OUTPUT:
        qSimilarity - NxN matrix containing correlation values between fibers
    """

    qDistances = _pairwiseQDistance_matrix(inputVTK, scalarData, scalarType, no_of_jobs)

    sigmasq = np.square(sigma)
    qSimilarity = distance.gausKernel_similarity(qDistances, sigmasq)

    qSimilarity = np.array(qSimilarity)

    return qSimilarity

def _degreeMatrix(inputMatrix):
    """ An internal function used to compute the Degree matrix, D """

    # Determine the degree matrix
    degMat = np.diag(np.sum(inputMatrix, 0))

    return degMat

def _cluster_to_rgb(data):

    """ Generate cluster color from first three components of data """

    colour = data[:, 0:3]

    # Normalize color
    colour_len = np.sqrt(np.sum(np.power(colour, 2), 1))
    colour = np.divide(colour.T, colour_len).T

    # Convert range from 0 to 255
    colour = 127.5 + (colour * 127.5)

    return colour.astype('int')

def _format_outputVTK(polyData, clusterIdx, colour, data):
    """ Output polydata with colours, cluster numbers and coordinates """

    dataColour = vtk.vtkUnsignedCharArray()
    dataColour.SetNumberOfComponents(3)
    dataColour.SetName('DataColour')

    clusterNumber = vtk.vtkIntArray()
    clusterNumber.SetName('ClusterNumber')

    for fidx in range(0, polyData.GetNumberOfLines()):
        dataColour.InsertNextTuple3(
                colour[clusterIdx[fidx], 0], colour[clusterIdx[fidx], 1], colour[clusterIdx[fidx], 2])
        clusterNumber.InsertNextTuple1(int(clusterIdx[fidx]))

    polyData.GetCellData().AddArray(dataColour)
    polyData.GetCellData().AddArray(clusterNumber)

    return polyData

def _weightedSimilarity(inputVTK, scalarDataList=[], scalarTypeList=[], scalarWeightList=[],
                                        sigma=1, no_of_jobs=1):
    """ Computes and returns a single weighted similarity matrix.
          Weight list should include weight for distance and sum to 1 """

    if ((scalarWeightList == []) & (scalarDataList != [])):
        print "\nNo weights given for provided measurements! Exiting..."
        return
    elif ((scalarDataList != []) & (scalarTypeList == [])):
        print "\nPlease also specify measurement(s) type. Exiting..."
        exit()
    elif (scalarDataList == []):
        print "\nNo measurements provided!"
        print "\nCalculating similarity based on geometry."
        wSimilarity = _pairwiseSimilarity_matrix(inputVTK, sigma, no_of_jobs)
    else:   # Calculate weighted similarity
        wSimilarity = _pairwiseSimilarity_matrix(inputVTK, sigma,
                                                                            no_of_jobs) * scalarWeightList[0]

        for i in range(len(scalarDataList)):
            similarity = _pairwiseQSimilarity_matrix(inputVTK, scalarDataList[i],
                scalarTypeList[i], sigma, no_of_jobs) * scalarWeightList[i+1]
            wSimilarity += similarity

        del similarity

    return wSimilarity
