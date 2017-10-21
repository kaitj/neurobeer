""" cluster.py

Module containing classes and functions used to cluster fibers and modify
parameters pertaining to clusters.

"""

import numpy as np
import os, scipy.cluster, sklearn.preprocessing
from joblib import Parallel, delayed
from joblib.pool import has_shareable_memory
from sys import exit

import fibers, distance, misc, prior
import vtk

def spectralClustering(inputVTK, scalarDataList=[], scalarTypeList=[], scalarWeightList=[],
                                    pts_per_fiber=20, k_clusters=200, sigma=0.4, saveAllSimilarity=False,
                                    saveWSimilarity=False, dirpath=None, verbose=0, no_of_jobs=1):
        """
        Clustering of fibers based on pairwise fiber similarity.
        See paper: "A tutorial on spectral clustering" (von Luxburg, 2007)

        If no scalar data provided, clustering performed based on geometry.
        First element of scalarWeightList should be weight placed for geometry, followed by order
        given in scalarTypeList These weights should sum to 1.0 (weights given as a decimal value).
        ex.  scalarDataList
              scalarTypeList = [FA, T1]
              scalarWeightList = [Geometry, FA, T1]

        INPUT:
            inputVTK - input polydata file
            scalarDataList - list containing scalar data for similarity measurements; defaults empty
            scalarTypeList - list containing scalar type for similarity measurements; defaults empty
            scalarWeightList - list containing scalar weights for similarity measurements; defaults empty
            pts_per_fiber - number of samples to take along each fiber
            k_clusters - number of clusters via k-means clustering; defaults 10 clusters
            sigma - width of Gaussian kernel; adjust to alter sensitivity; defaults 0.4
            saveAllSimilarity - flag to save all individual similarity matrices computed; defaults False
            saveWSimilarity - flag to save weighted similarity matrix computed; defaults False
            dirpath - directory to store files; defaults None
            verbose - verbosity of function; defaults 0
            no_of_jobs - cores to use to perform computation; defaults 1

        OUTPUT:
            outputPolydata - polydata containing information from clustering to be written into VTK
            clusterIdx - array containing cluster that each fiber belongs to
            fiberData - tree containing spatial and quantitative information of fibers
        """
        if dirpath is None:
            dirpath = os.getcwd()

        matpath = dirpath + '/matrices'
        if not os.path.exists(matpath):
            os.makedirs(matpath)

        noFibers = inputVTK.GetNumberOfLines()
        if noFibers == 0:
            print "\nERROR: Input data has 0 fibers!"
            raise ValueError
        elif verbose == 1:
            print "\nStarting clustering..."
            print "No. of fibers:", noFibers
            print "No. of clusters:", k_clusters

        fiberData = fibers.FiberTree()
        fiberData.convertFromVTK(inputVTK, pts_per_fiber, verbose)
        for i in range(len(scalarTypeList)):
            fiberData.addScalar(inputVTK, scalarDataList[i], scalarTypeList[i], pts_per_fiber)

        # 1. Compute similarty matrix
        W = _pairwiseWeightedSimilarity(fiberData, scalarTypeList, scalarWeightList,
                                                    sigma, saveAllSimilarity, pts_per_fiber, matpath, no_of_jobs)

        if saveWSimilarity is True:
            misc.saveMatrix(matpath, W, 'Weighted')

        # 2. Compute degree matrix
        W = np.dot(W, W.T)  # Computes correlation matrix from similarity
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
        misc.saveEig(dirpath, eigval, eigvec)

        Dscale = np.divide(1, np.sqrt(np.sum(D, 0)))
        em_vec = np.multiply(Dscale, eigvec)

        # 6. Compute information for clustering using "N" number of smallest eigenvalues
        # Skips first eigenvector, no information obtained
        if k_clusters > eigvec.shape[0]:
            print '\nNumber of user selected clusters greater than number of eigenvectors.'
            print 'Clustering with maximum number of available eigenvectors.'
            U = em_vec[:, 1:em_vec.shape[0]]
        elif k_clusters == eigvec.shape[0]:
            U = em_vec[:, 1:k_clusters]
        else:
            U = em_vec[:, 1:k_clusters + 1]
        U = U.real

        # 7. Find clusters using K-means clustering
        centroids, clusterIdx = scipy.cluster.vq.kmeans2(U, k_clusters, iter=20, minit='points')
        centroids, clusterIdx = _sortLabel(centroids, clusterIdx)
        fiberData.addClusterInfo(clusterIdx, centroids)

        if k_clusters <= 1:
            print "\nNot enough eigenvectors selected!"
            raise ValueError
        elif k_clusters == 2:
            temp = eigvec[:, 0:3]
            temp = temp.astype('float')
            colour = _cluster_to_rgb(temp)
            del temp
        else:
            colour = _cluster_to_rgb(centroids)

        # 8. Return results
        # Create model with user / default number of chosen samples along fiber
        outputData = fiberData.convertToVTK()

        outputPolydata = _format_outputVTK(outputData, clusterIdx, colour, centroids)

        # 9. Also add measurements from those used to cluster
        for i in range(len(scalarTypeList)):
            outputPolydata = addScalarToVTK(outputPolydata, fiberData, scalarTypeList[i])

        return outputPolydata, clusterIdx, fiberData

def spectralPriorCluster(inputVTK, priorVTK, scalarDataList=[], scalarTypeList=[],
                                    scalarWeightList=[], sigma=0.4, saveAllSimilarity=False,
                                    saveWSimilarity=False, dirpath=None, verbose=0, no_of_jobs=1):
        """
        Clustering of fibers based on pairwise fiber similarity using previously clustered fibers.
        See paper: "A tutorial on spectral clustering" (von Luxburg, 2007)

        If no scalar data provided, clustering performed based on geometry.
        First element of scalarWeightList should be weight placed for geometry, followed by order
        given in scalarTypeList These weights should sum to 1.0 (weights given as a decimal value).
        ex. scalarDataList
              scalarTypeList = [FA, T1]
              scalarWeightList = [Geometry, FA, T1]

        INPUT:
            inputVTK - input polydata file
            priorVTK - prior polydata file
            scalarDataList - list containing scalar data for similarity measurements; defaults empty
            scalarTypeList - list containing scalar type for similarity measurements; defaults empty
            scalarWeightList - list containing scalar weights for similarity measurements; defaults empty
            sigma - width of Gaussian kernel; adjust to alter sensitivity; defaults 0.4
            saveAllSimilarity - flag to save all individual similarity matrices computed; defaults False
            saveWSimilarity - flag to save weighted similarity matrix computed; defaults False
            dirpath - directory to store files; defaults None
            verbose - verbosity of function; defaults 0
            no_of_jobs - cores to use to perform computation; defaults 1

        OUTPUT:
            outputPolydata - polydata containing information from clustering to be written into VTK
            clusterIdx - array containing cluster that each fiber belongs to
            fiberData - tree containing spatial and quantitative information of fibers
        """

        priorData, priorCentroids = prior.load(priorVTK)

        k_clusters = len(priorCentroids)
        pts_per_fiber = int(priorData.pts_per_fiber)

        noFibers = inputVTK.GetNumberOfLines()
        nopriorFibers = int(priorData.no_of_fibers)
        if noFibers == 0 or nopriorFibers == 0:
            print "\nERROR: Input data(s) has 0 fibers!"
            raise ValueError
        elif verbose == 1:
            print "\nStarting clustering..."
            print "No. of fibers:", noFibers
            print "No. of clusters:", k_clusters

        fiberData = fibers.FiberTree()
        fiberData.convertFromVTK(inputVTK, pts_per_fiber, verbose)
        for i in range(len(scalarTypeList)):
            fiberData.addScalar(inputVTK, scalarDataList[i], scalarTypeList[i], pts_per_fiber)

        # 1. Compute similarty matrix
        W = _priorWeightedSimilarity(fiberData, priorData, scalarTypeList, scalarWeightList,
                                                    sigma, saveAllSimilarity, pts_per_fiber, dirpath, no_of_jobs)
        V = np.dot(W, W.T)  # Computes low order approximation from left SVD

        if dirpath is None:
            dirpath = os.getcwd()
            if not os.path.exists(dirpath + '/eigval.npy') or os.path.exists(dirpath + '/eigvec.npy'):
                print "\nMissing eigenvalue or eigenvector binary file"
                raise "I/O Error"

            eigval, eigvec = prior.loadEig(dirpath, 'eigval.npy', 'eigvec.npy')

        if saveWSimilarity is True:
            matpath = dirpath + '/matrices'
            if not os.path.exists(matpath):
                os.makedirs(matpath)

            misc.saveMatrix(matpath, W, 'Weighted')
            misc.saveMatrix(matpath, V, 'Approximation')

        # 2. Compute degree matrix
        D = _degreeMatrix(V)
        Dscale = np.divide(1, np.sqrt(np.sum(D, 0)))

        # 3. Compute embedding vector
        em_vec = np.multiply(Dscale, eigvec)

        # 6. Compute information for clustering using "N" number of smallest eigenvalues
        # Skips first eigenvector, no information obtained
        if k_clusters > eigvec.shape[0]:
            print 'Number of user selected clusters greater than number of eigenvectors.'
            print 'Clustering with maximum number of available eigenvectors.'
            U = em_vec[:, 1:em_vec.shape[0]]
        elif k_clusters == eigvec.shape[0]:
            U = em_vec[:, 1:k_clusters]
        else:
            U = em_vec[:, 1:k_clusters + 1]
        U = U.real

        # 7. Find clusters using K-means clustering
        clusterIdx, dist = scipy.cluster.vq.vq(U, priorCentroids)
        fiberData.addClusterInfo(clusterIdx, priorCentroids)

        if k_clusters <= 1:
            print('Not enough eigenvectors selected!')
            raise ValueError
        elif k_clusters == 2:
            temp = eigvec[:, 0:3]
            temp = temp.astype('float')
            colour = _cluster_to_rgb(temp)
            del temp
        else:
            colour = _cluster_to_rgb(priorCentroids)

        # 8. Return results
        # Create model with user / default number of chosen samples along fiber
        outputData = fiberData.convertToVTK()

        outputPolydata = _format_outputVTK(outputData, clusterIdx, colour, priorCentroids)

        # 9. Also add measurements from those used to cluster
        for i in range(len(scalarTypeList)):
            outputPolydata = addScalarToVTK(outputPolydata, fiberData, scalarTypeList[i])

        return outputPolydata, clusterIdx, fiberData

def addScalarToVTK(polyData, fiberTree, scalarType, fidxes=None):
    """
    Add scalar to all polydata points to be converted to .vtk file.
    This function is different from scalars.addScalar, which only considers point
    used in sampling of fiber.

    INPUT:
        polyData - polydata for scalar measurements to be added to
        fiberTree - the tree containing polydata information
        scalarType - type of quantitative measurement to be aded to polydata
        fidxes - array with fiber indices pertaining to scalar data of extracted fibers; default none

    OUTPUT:
        polydata - updated polydata with quantitative information
    """

    data = vtk.vtkFloatArray()
    data.SetName(scalarType.split('/', -1)[-1])

    if fidxes is None:
        for fidx in range(0, polyData.GetNumberOfLines()):
            for pidx in range(0, fiberTree.pts_per_fiber):
                scalarValue = fiberTree.fiberTree[fidx][pidx][scalarType]
                data.InsertNextValue(float(scalarValue))
    else:
        for fidx in fidxes:
            for pidx in range(0, fiberTree.pts_per_fiber):
                scalarValue = fiberTree.fiberTree[fidx][pidx][scalarType]
                data.InsertNextValue(float(scalarValue))

    polyData.GetPointData().AddArray(data)

    return polyData

def extractCluster(inputVTK, clusterIdx, label, pts_per_fiber):
    """
    Extracts a cluster corresponding to the label provided.

    INPUT:
        inputVTK - polydata to extract cluster from
        clusterIdx - labels pertaining to fibers of inputVTK
        label - label of cluster to be extracted
        pts_per_fiber - number of samples to take along fiber

    OUTPUT:
        polyData - extracted cluster in polydata format; no information is retained
    """
    fiberTree = fibers.FiberTree()
    fiberTree.convertFromVTK(inputVTK, pts_per_fiber)

    cluster = fiberTree.getFibers(np.where(clusterIdx == label)[0])
    cluster = fibers.convertFromTuple(cluster)
    polyData = cluster.convertToVTK()

    return polyData

def _pairwiseDistance_matrix(fiberTree, pts_per_fiber, no_of_jobs):
    """ *INTERNAL FUNCTION*
    Used to compute an NxN distance matrix for all fibers (N) in the input data.

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        pts_per_fiber - number of samples along a fiber
        no_of_jobs - cores to use to perform computation

    OUTPUT:
        distances - NxN matrix containing distances between fibers
    """

    distances = Parallel(n_jobs=no_of_jobs, backend="threading")(
            delayed(distance.fiberDistance, has_shareable_memory)(fiberTree.getFiber(fidx),
                    fiberTree.getFibers(range(fiberTree.no_of_fibers)))
            for fidx in range(0, fiberTree.no_of_fibers))

    distances = np.array(distances)
    # Normalize between 0 and 1
    distances = sklearn.preprocessing.MinMaxScaler().fit_transform(distances)

    if np.diag(distances).all() != 0.0:
        print('Diagonals in distance matrix are not equal to 0')
        exit()

    return distances

def _pairwiseSimilarity_matrix(fiberTree, sigma, pts_per_fiber, no_of_jobs):
    """ *INTERNAL FUNCTION*
    Computes an NxN similarity matrix for all fibers (N) in the input data.

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        sigma - width of Gaussian kernel; adjust to alter
        pts_per_fiber - number of samples along a fiber
        no_of_jobs - cores to use to perform computation

    OUTPUT:
        similarity - NxN matrix containing similarity between fibers based on geometry
    """

    distances = _pairwiseDistance_matrix(fiberTree, pts_per_fiber, no_of_jobs)

    sigmasq = np.square(sigma)
    similarities = distance.gausKernel_similarity(distances, sigmasq)

    similarities = np.array(similarities)

    if np.diag(similarities).all() != 1.0:
        print('Diagonals in similarity matrix are not equal to 1')
        exit()

    return similarities

def _pairwiseQDistance_matrix(fiberTree, scalarType, pts_per_fiber, no_of_jobs):
    """ *INTERNAL FUNCTION*
    Computes the "pairwise distance" between quantitative points along a fiber.

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        scalarType - type of quantitative measurements to be used for computation
        pts_per_fiber - number of sample along a fiber
        no_of_jobs - cores to use to perform computation

    OUTPUT:
        qDistances - NxN matrix containing pairwise distances between fibers
    """

    no_of_fibers = fiberTree.no_of_fibers

    qDistances = Parallel(n_jobs=no_of_jobs, backend="threading")(
            delayed(distance.scalarDistance, has_shareable_memory)(
                fiberTree.getScalar(fidx, scalarType),
                fiberTree.getScalars(range(no_of_fibers), scalarType))
            for fidx in range(0, no_of_fibers)
    )

    qDistances = np.array(qDistances)

    # Normalize distance measurements
    qDistances = sklearn.preprocessing.MinMaxScaler().fit_transform(qDistances)

    if np.diag(qDistances).all() != 0.0:
        print "Diagonals in distance matrix are not equal to 0"
        exit()

    return qDistances

def _pairwiseQSimilarity_matrix(fiberTree, scalarType, sigma, pts_per_fiber,
                                                      no_of_jobs):
    """ *INTERNAL FUNCTION*
    Computes the similarity between quantitative points along a fiber.

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        scalarType - type of quantitative measurements to be used for computation
        sigma - width of Gaussian kernel; adjust to alter sensitivity
        no_of_jobs - cores to use to perform computation
        pts_per_fiber - number of samples along a fiber

    OUTPUT:
        qSimilarity - NxN matrix containing similarity of quantitative measurements between fibers
    """

    qDistances = _pairwiseQDistance_matrix(fiberTree, scalarType, pts_per_fiber, no_of_jobs)

    sigmasq = np.square(sigma)
    qSimilarity = distance.gausKernel_similarity(qDistances, sigmasq)

    qSimilarity = np.array(qSimilarity)

    if np.diag(qSimilarity).all() != 1.0:
        print "Diagonals in similarity marix are not equal to 1"
        exit()

    return qSimilarity

def _priorDistance_matrix(fiberTree, priorTree, pts_per_fiber, no_of_jobs):
    """ *INTERNAL FUNCTION*
    Used to compute an distance matrix for all fibers (N) in the input data through
    comparison with previously clustered data

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        priorTree - tree containing spatial and quantitative info from prev. clustered fibers
        pts_per_fiber - number of samples along a fiber
        no_of_jobs - cores to use to perform computation

    OUTPUT:
        distances - matrix containing distances between fibers
    """

    distances = Parallel(n_jobs=no_of_jobs, backend="threading")(
            delayed(distance.fiberDistance, has_shareable_memory)(fiberTree.getFiber(fidx),
                    priorTree.getFibers(range(priorTree.no_of_fibers)))
            for fidx in range(0, fiberTree.no_of_fibers))

    distances = np.array(distances)

    # Normalize between 0 and 1
    distances = sklearn.preprocessing.MinMaxScaler().fit_transform(distances)

    return distances

def _priorSimilarity_matrix(fiberTree, priorTree, sigma, pts_per_fiber, no_of_jobs):
    """ *INTERNAL FUNCTION*
    Computes a similarity matrix for all fibers (N) in the input data to previously clustered fibers

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        priorTree - tree containing spatial and quantitative info from previously clustered fibers
        sigma - width of Gaussian kernel; adjust to alter
        pts_per_fiber - number of samples along a fiber
        no_of_jobs - cores to use to perform computation

    OUTPUT:
        similarity - matrix containing similarity between fibers based on geometry
    """

    distances = _priorDistance_matrix(fiberTree, priorTree, pts_per_fiber, no_of_jobs)

    sigmasq = np.square(sigma)
    similarities = distance.gausKernel_similarity(distances, sigmasq)

    similarities = np.array(similarities)

    return similarities

def _priorQDistance_matrix(fiberTree, priorTree, scalarType, pts_per_fiber, no_of_jobs):
    """ *INTERNAL FUNCTION*
    Computes the "pairwise distance" between quantitative points along a fiber and previously
    clustered fibers

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        priorTree - tree containing information on previously clustered fibers
        scalarType - type of quantitative measurements to be used for computation
        pts_per_fiber - number of sample along a fiber
        no_of_jobs - cores to use to perform computation

    OUTPUT:
        qDistances - matrix containing pairwise distances between fibers
    """

    qDistances = Parallel(n_jobs=no_of_jobs, backend="threading")(
            delayed(distance.scalarDistance, has_shareable_memory)(
                fiberTree.getScalar(fidx, scalarType),
                priorTree.getScalars(range(priorTree.no_of_fibers), scalarType))
            for fidx in range(0, fiberTree.no_of_fibers)
    )

    qDistances = np.array(qDistances)

    # Normalize distance measurements
    qDistances = sklearn.preprocessing.MinMaxScaler().fit_transform(qDistances)

    return qDistances

def _priorQSimilarity_matrix(fiberTree, priorTree, scalarType, sigma, pts_per_fiber,
                                                      no_of_jobs):
    """ *INTERNAL FUNCTION*
    Computes the similarity between quantitative points along a fiber and previously clustered
    fibers

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        priorsTree - tree containing information on previously clustered fibers
        scalarType - type of quantitative measurements to be used for computation
        sigma - width of Gaussian kernel; adjust to alter sensitivity
        no_of_jobs - cores to use to perform computation
        pts_per_fiber - number of samples along a fiber

    OUTPUT:
        qSimilarity - matrix containing similarity of quantitative measurements between fibers
    """

    qDistances = _priorQDistance_matrix(fiberTree, priorTree, scalarType, pts_per_fiber, no_of_jobs)

    sigmasq = np.square(sigma)
    qSimilarity = distance.gausKernel_similarity(qDistances, sigmasq)

    qSimilarity = np.array(qSimilarity)

    if np.diag(qSimilarity).all() != 1.0:
        print "Diagonals in similarity marix are not equal to 1"
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
    colourMag = np.sqrt(np.sum(np.square(colour), 1))
    colour = np.divide(colour.T, colourMag).T

    # Convert range from 0 to 255
    colour = 127.5 + (colour * 127.5)

    return colour.astype('int')

def _format_outputVTK(polyData, clusterIdx, colour, centroids):
    """ *INTERNAL FUNCTION*
    Formats polydata with cluster index and colour.

    INPUT:
        polyData - polydata for information to be applied to
        clusterIdx - cluster indices to be applied to each fiber within the polydata model
        colour - colours to be applied to each fiber within the polydata model
        centroid - centroid location to associated with each cluster

    OUTPUT:
        polyData - updated polydata with cluster and colour information
    """

    dataColour = vtk.vtkUnsignedCharArray()
    dataColour.SetNumberOfComponents(3)
    dataColour.SetName('Colour')

    clusterLabel = vtk.vtkIntArray()
    clusterLabel.SetNumberOfComponents(1)
    clusterLabel.SetName('ClusterLabel')

    centroid = vtk.vtkFloatArray()
    centroid.SetNumberOfComponents(centroids.shape[1])
    centroid.SetName('Centroid')

    for fidx in range(0, polyData.GetNumberOfLines()):
        dataColour.InsertNextTuple3(
                colour[clusterIdx[fidx], 0], colour[clusterIdx[fidx], 1], colour[clusterIdx[fidx], 2])
        clusterLabel.InsertNextTuple1(int(clusterIdx[fidx]))
        centroid.InsertNextTuple(centroids[clusterIdx[fidx], :])

    polyData.GetCellData().AddArray(dataColour)
    polyData.GetCellData().AddArray(clusterLabel)
    polyData.GetCellData().AddArray(centroid)

    return polyData

def _pairwiseWeightedSimilarity(fiberTree, scalarTypeList=[], scalarWeightList=[],
                                        sigma=0.4, saveAllSimilarity=False, pts_per_fiber=20, dirpath=None,
                                        no_of_jobs=1):
    """ *INTERNAL FUNCTION*
    Computes and returns a single weighted similarity matrix.
    Weight list should include weight for distance and sum to 1.

    INPUT:
        fiberTree - tree containing scalar data for similarity measurements
        scalarTypeList - list containing scalar type for similarity measurements; defaults empty
        scalarWeightList - list containing scalar weights for similarity measurements; defaults empty
        sigma - width of Gaussian kernel; adjust to alter sensitivity; defaults 0.4
        saveAllSimilarity - flag to save all individual similarity matrices computed; defaults 0 (off)
        dirpath - directory to store similarity matrices
        no_of_jobs - cores to use to perform computation; defaults 1

    OUTPUT:
        wSimilarity - matrix containing the computed weighted similarity
    """

    if ((scalarWeightList == []) & (scalarTypeList != [])):
        print "\nNo weights given for provided measurements! Exiting..."
        exit()

    elif ((scalarWeightList != []) & (scalarTypeList == [])):
        print "\nPlease also specify measurement(s) type. Exiting..."
        exit()

    elif ((scalarWeightList == [])) & ((scalarTypeList == [])):
        print "\nNo measurements provided!"
        print "\nCalculating similarity based on geometry."
        wSimilarity = _pairwiseSimilarity_matrix(fiberTree, sigma, pts_per_fiber, no_of_jobs)

        if dirpath is None:
            dirpath = os.getcwd()

        misc.saveMatrix(dirpath, wSimilarity, 'Geometry')

    else:   # Calculate weighted similarity

        if np.sum(scalarWeightList) != 1.0:
            print '\nWeights given do not sum 1. Exiting...'
            exit()

        wSimilarity = _pairwiseSimilarity_matrix(fiberTree, sigma, pts_per_fiber,
                                                                            no_of_jobs)

        if saveAllSimilarity is True:
            if dirpath is None:
                dirpath = os.getcwd()

            matrixType = scalarTypeList[0].split('/', -1)[-1]
            matrixType = matrixType[:-2] + 'geometry'
            misc.saveMatrix(dirpath, wSimilarity, matrixType)

        wSimilarity = wSimilarity * scalarWeightList[0]

        for i in range(len(scalarTypeList)):
            similarity = _pairwiseQSimilarity_matrix(fiberTree,
                scalarTypeList[i], sigma, pts_per_fiber, no_of_jobs)

            if saveAllSimilarity is True:
                if dirpath is None:
                    dirpath = os.getcwd()

                matrixType = scalarTypeList[i].split('/', -1)[-1]
                misc.saveMatrix(dirpath, similarity, matrixType)

            wSimilarity += similarity * scalarWeightList[i+1]

        del similarity

    if np.diag(wSimilarity).all() != 1.0:
        print "Diagonals of weighted similarity are not equal to 1"
        exit()

    return wSimilarity

def _priorWeightedSimilarity(fiberTree, priorTree, scalarTypeList=[], scalarWeightList=[],
                                        sigma=0.4, saveAllSimilarity=False, pts_per_fiber=20, dirpath=None,
                                        no_of_jobs=1):
    """ *INTERNAL FUNCTION*
    Computes and returns a single weighted similarity matrix.
    Weight list should include weight for distance and sum to 1.

    INPUT:
        fiberTree - tree containing scalar data for similarity measurements
        priorTree - tree containing previously clustered tract information
        scalarTypeList - list containing scalar type for similarity measurements; defaults empty
        scalarWeightList - list containing scalar weights for similarity measurements; defaults empty
        sigma - width of Gaussian kernel; adjust to alter sensitivity; defaults 0.4
        saveAllSimilarity - flag to save all individual similarity matrices computed; defaults 0 (off)
        dirpath - directory to store similarity matrices
        no_of_jobs - cores to use to perform computation; defaults 1

    OUTPUT:
        wSimilarity - matrix containing the computed weighted similarity
    """

    if ((scalarWeightList == []) & (scalarTypeList != [])):
        print "\nNo weights given for provided measurements! Exiting..."
        exit()

    elif ((scalarWeightList != []) & (scalarTypeList == [])):
        print "\nPlease also specify measurement(s) type. Exiting..."
        exit()

    elif ((scalarWeightList == [])) & ((scalarTypeList == [])):
        print "\nNo measurements provided!"
        print "\nCalculating similarity based on geometry."
        wSimilarity = _priorSimilarity_matrix(fiberTree, priorTree, sigma, pts_per_fiber, no_of_jobs)

        if dirpath is None:
            dirpath = os.getcwd()
        else:
            if not os.path.exists(dirpath):
                os.makedirs(dirpath)

        misc.saveMatrix(dirpath, wSimilarity, 'Geometry')

    else:   # Calculate weighted similarity

        if np.sum(scalarWeightList) != 1.0:
            print '\nWeights given do not sum 1. Exiting...'
            exit()

        wSimilarity = _priorSimilarity_matrix(fiberTree, priorTree, sigma, pts_per_fiber,
                                                                            no_of_jobs)

        if saveAllSimilarity is True:
            if dirpath is None:
                dirpath = os.getcwd()

            matrixType = scalarTypeList[0].split('/', -1)[-1]
            matrixType = matrixType[:-2] + 'geometry'
            misc.saveMatrix(dirpath, wSimilarity, matrixType)

        wSimilarity = wSimilarity * scalarWeightList[0]

        for i in range(len(scalarTypeList)):
            similarity = _priorQSimilarity_matrix(fiberTree, priorTree,
                scalarTypeList[i], sigma, pts_per_fiber, no_of_jobs)

            if saveAllSimilarity is True:
                if dirpath is None:
                    dirpath = os.getcwd()

                matrixType = scalarTypeList[i].split('/', -1)[-1]
                misc.saveMatrix(dirpath, similarity, matrixType)

            wSimilarity += similarity * scalarWeightList[i+1]

        del similarity

    if np.diag(wSimilarity).all() != 1.0:
        print "\nDiagonals of weighted similarity are not equal to 1"
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
        newCentroids[i, :] = centroids[sortedClusters[i]]

    return newCentroids, newClusters
