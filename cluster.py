""" cluster.py

Module containing classes and functions used to cluster fibers and modify
parameters pertaining to clusters.

"""

import numpy as np
import os, scipy.cluster
from joblib import Parallel, delayed
from sys import exit

import fibers, distance, scalars, misc
import vtk

import sklearn.cluster, sklearn.preprocessing

def spectralClustering(inputVTK, scalarDataList=[], scalarTypeList=[], scalarWeightList=[],
                                    k_clusters=10, sigma=0.4, saveAllSimilarity=0,
                                    saveWSimilarity=0, no_of_jobs=1):
        """
        Clustering of fibers based on pairwise fiber similarity.
        See paper: "A tutorial on spectral clustering" (von Luxburg, 2007)

        If no scalar lists are provided, clustering performed using tract geometry.
        First element of scalarWeightList should be weight placed for geometry, followed by order
        given in scalarTypeList These weights should sum to 1.0 (weights given as a decimal value).
        ex. scalarDataList = [FA, T1]
              scalarTypeList = [FA, T1]
              scalarWeightList = [Geometry, FA, T1]

        INPUT:
            inputVTK - input polydata file
            scalarDataList - list containing scalar data for similarity measurements; defaults empty
            scalarTypeList - list containing scalar type for similarity measurements; defaults empty
            scalarWeightList - list containing scalar weights for similarity measurements; defaults empty
            k_clusters - number of clusters via k-means clustering; defaults 10 clusters
            sigma - width of Gaussian kernel; adjust to alter sensitivity; defaults 0.4
            saveAllSimilarity - flag to save all individual similarity matrices computed; defaults 0 (off)
            saveWSimilarity - flag to save weighted similarity matrix computed; defaults 0 (off)
            no_of_jobs - cores to use to perform computation; defaults 1

        OUTPUT: outputPolydata, clusterIdx, colour, centroids
            outputPolydata - polydata containing information from clustering to be written into VTK
            clusterIdx - array containing cluster that each fiber belongs to
            colour - array containing the RGB value assigned to the scalar
            centroids - array containing the centroids for each cluster

        NOTE: For most consistent results, k_clusters => no_of_eigvec
        """

        no_of_eigvec = k_clusters

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
                                                    sigma, saveAllSimilarity, no_of_jobs)

        if saveWSimilarity == 1:
            # Temporary solution, need to decide how to store location to save matrix
            dirpath = raw_input("Enter path to save weighted similarity matrix: ")
            if not os.path.exists(dirpath):
                os.makedirs(dirpath)

            misc.saveMatrix(dirpath, W, 'Weighted')

        # 2. Compute degree matrix
        D = _degreeMatrix(W)

        # 3. Compute unnormalized Laplacian
        L = D - W

        # 4. Compute normalized Laplacian (random-walk)
        Lrw = np.dot(np.diag(np.divide(1, np.sum(D, 0))), L)

        # 5. Compute eigenvalues and eigenvectors of generalized eigenproblem
        # Sort by ascending eigenvalue
        eigval, eigvec = np.linalg.eig(Lrw)
        idx = eigval.argsort()
        eigval, eigvec = eigval[idx], eigvec[:, idx]

        # 6. Compute information for clustering using "N" number of smallest eigenvalues
        U = eigvec[:, 0:no_of_eigvec]
        U = U.astype('float')

        # 7. Find clusters using K-means clustering
        centroids, clusterIdx = scipy.cluster.vq.kmeans2(U, k_clusters, iter=20, minit='points')
        centroids, clusterIdx = _sortLabel(centroids, clusterIdx)

        if no_of_eigvec <= 1:
            print('Not enough eigenvectors selected!')
        elif no_of_eigvec == 2:
            temp = eigvec[:, 0:3]
            temp = temp.astype('float')
            colour = _cluster_to_rgb(temp)
            del temp
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
    """
    Add scalar to polydata points to be converted to .vtk file.
    This function is different from scalars.addScalar, which only considers point
    used in sampling of fiber.

    INPUT:
        polyData - polydata for scalar measurements to be added to
        scalarData - the array containing the quantitative measurement  to be added to the polydata
        scalarType - type of quantitative measurement to be aded to polydata

    OUTPUT:
        polydata - updated polydata with quantitative information
    """

    data = vtk.vtkFloatArray()
    data.SetName(scalarType.split('/', -1)[-1])

    for pidx in range(0, polyData.GetNumberOfPoints()):
        data.InsertNextValue(float(scalarData[pidx]))

    polyData.GetPointData().AddArray(data)

    return polyData

def _pairwiseDistance_matrix(inputVTK, no_of_jobs):
    """ *INTERNAL FUNCTION*
    Used to compute an NxN distance matrix for all fibers (N) in the input data.

    INPUT:
        inputVTK - input polydata file
        no_of_jobs - cores to use to perform computation

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

    if np.diag(distances).all() != 0.0:
        print('Diagonals in distance matrix are not equal to 0')
        exit()

    return distances

def _pairwiseSimilarity_matrix(inputVTK, sigma, no_of_jobs):
    """ *INTERNAL FUNCTION*
    Computes an NxN similarity matrix for all fibers (N) in the input data.

    INPUT:
        inputVTK - input polydata file
        sigma - width of Gaussian kernel; adjust to alter sensitivity
        no_of_jobs - cores to use to perform computation

    OUTPUT:
        similarity - NxN matrix containing similarity between fibers based on geometry
    """

    distances = _pairwiseDistance_matrix(inputVTK, no_of_jobs)

    sigmasq = np.square(sigma)
    similarities = distance.gausKernel_similarity(distances, sigmasq)

    similarities = np.array(similarities)

    if np.diag(similarities).all() != 1.0:
        print('Diagonals in similarity matrix are not equal to 1')
        exit()

    return similarities

def _pairwiseQDistance_matrix(inputVTK, scalarData, scalarType, no_of_jobs):
    """ *INTERNAL FUNCTION*
    Computes the "pairwise distance" between quantitative points along a fiber.

    INPUT:
        inputVTK - input polydata file
        scalarData - array containing quantitative measurements to be used for computation
        scalarType - type of quantitative measurements to be used for computation
        no_of_jobs - cores to use to perform computation

    OUTPUT:
        qDistances - NxN matrix containing pairwise distances between fibers
    """

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

    if np.diag(qDistances).all() != 0.0:
        print('Diagonals in distance matrix are not equal to 0')
        exit()

    return qDistances

def _pairwiseQSimilarity_matrix(inputVTK, scalarData, scalarType, sigma, no_of_jobs):
    """ *INTERNAL FUNCTION*
    Computes the similarity between quantitative points along a fiber.

    INPUT:
        inputVTK - input polydata file
        scalarData - array containing quantitative measurements to be used for computation
        scalarType - type of quantitative measurements to be used for computation
        sigma - width of Gaussian kernel; adjust to alter sensitivity
        no_of_jobs - cores to use to perform computation

    OUTPUT:
        qSimilarity - NxN matrix containing similarity of quantitative measurements between fibers
    """

    qDistances = _pairwiseQDistance_matrix(inputVTK, scalarData, scalarType, no_of_jobs)

    sigmasq = np.square(sigma)
    qSimilarity = distance.gausKernel_similarity(qDistances, sigmasq)

    qSimilarity = np.array(qSimilarity)

    if np.diag(qSimilarity).all() != 1.0:
        print('Diagonals in similarity marix are not equal to 1')
        exit()

    return qSimilarity

def _degreeMatrix(inputMatrix):
    """ *INTERNAL FUNCTION*
    Computes the degree matrix, D.

    INPUT:
        inputMatrix - adjacency matrix to be used for computation

    OUTPUT:
        degMat - degree matrix to be used to compute Laplacian matrix
    """

    # Determine the degree matrix
    degMat = np.diag(np.sum(inputMatrix, 0))

    return degMat

def _cluster_to_rgb(data):
    """ *INTERNAL FUNCTION*
    Generate cluster color from first three components of data

    INPUT:
        data - information used to calculate RGB colours; typically eigenvectors or centroids are used

    OUTPUT:
        colour - array containing the RGB values to colour clusters
    """

    colour = data[:, 0:3]

    # Normalize color
    colour_len = np.sqrt(np.sum(np.power(colour, 2), 1))
    colour = np.divide(colour.T, colour_len).T

    # Convert range from 0 to 255
    colour = 127.5 + (colour * 127.5)

    return colour.astype('int')

def _format_outputVTK(polyData, clusterIdx, colour):
    """ *INTERNAL FUNCTION*
    Formats polydata with cluster index and colour.

    INPUT:
        polyData - polydata for information to be applied to
        clusterIdx - cluster indices to be applied to each fiber within the polydata model
        colour - colours to be applied to each fiber within the polydata model

    OUTPUT:
        polyData - updated polydata with cluster and colour information
    """

    dataColour = vtk.vtkUnsignedCharArray()
    dataColour.SetNumberOfComponents(3)
    dataColour.SetName('Colour')

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
                                        sigma=0.4, saveAllSimilarity=0, no_of_jobs=1):
    """ *INTERNAL FUNCTION*
    Computes and returns a single weighted similarity matrix.
    Weight list should include weight for distance and sum to 1.

    INPUT:
        inputVTK - input polydata file
        scalarDataList - list containing scalar data for similarity measurements; defaults empty
        scalarTypeList - list containing scalar type for similarity measurements; defaults empty
        scalarWeightList - list containing scalar weights for similarity measurements; defaults empty
        sigma - width of Gaussian kernel; adjust to alter sensitivity; defaults 0.4
        saveAllSimilarity - flag to save all individual similarity matrices computed; defaults 0 (off)
        no_of_jobs - cores to use to perform computation; defaults 1

    OUTPUT:
        wSimilarity - matrix containing the computed weighted similarity
    """

    if ((scalarWeightList == []) & (scalarDataList != [])):
        print "\nNo weights given for provided measurements! Exiting..."
        exit()

    elif ((scalarDataList != []) & (scalarTypeList == [])):
        print "\nPlease also specify measurement(s) type. Exiting..."
        exit()

    elif (scalarDataList == []):
        print "\nNo measurements provided!"
        print "\nCalculating similarity based on geometry."
        wSimilarity = _pairwiseSimilarity_matrix(inputVTK, sigma, no_of_jobs)

        if (saveAllSimilarity == 1):
            dirpath = raw_input("Enter path to store similarity matrix: ")
            if not os.path.exists(dirpath):
                os.makedirs(dirpath)

            misc.saveMatrix(dirpath, wSimilarity, 'Geometry')

    else:   # Calculate weighted similarity

        checksum = np.sum(scalarWeightList)
        if checksum > 1.0:
            print '\nWeights given sum to greater than 1. Exiting...'
            exit()

        wSimilarity = _pairwiseSimilarity_matrix(inputVTK, sigma,
                                                                            no_of_jobs) * scalarWeightList[0]

        for i in range(len(scalarDataList)):
            similarity = _pairwiseQSimilarity_matrix(inputVTK, scalarDataList[i],
                scalarTypeList[i], sigma, no_of_jobs)

            if saveAllSimilarity == 1:
                dirpath = raw_input("Enter path to store similarity matrix: ")
                if not os.path.exists(dirpath):
                    os.makedirs(dirpath)

                matrixType = scalarTypeList[i].split('/', -1)[-1]
                misc.saveMatrix(dirpath, similarity, matrixType, 'weighted')

            wSimilarity += similarity * scalarWeightList[i+1]

        del similarity

    if np.diag(wSimilarity).all() != 1.0:
        print('Diagonals of weighted similarity are not equal to 1')
        exit()

    return wSimilarity

def _sortLabel(centroids, clusterIdx):
    """ *INTERNAL FUNCTION*
    Sort the cluster label by fiber count.

    INPUT:
        centroids - array of centroids to be sorted
        clusterIdx - array containing cluster indices to be sorted

    OUTPUT:
        newCentroids - array of sorted centroids
        newClusters - array of sorted clusters
    """

    uniqueClusters, countClusters = np.unique(clusterIdx, return_counts=True)
    sortedClusters = np.argsort(-countClusters)

    newClusters = np.copy(clusterIdx)
    newCentroids = np.copy(centroids)

    for i in range(len(sortedClusters)):
        newIdx = np.where(sortedClusters == i)
        newClusters[clusterIdx == i] = newIdx[0][0]
        newCentroids[i, :] = centroids[newIdx[0][0], :]

    return newCentroids, newClusters
