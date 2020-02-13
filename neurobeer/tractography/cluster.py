""" cluster.py

Module containing classes and functions used to cluster fibers and modify
parameters pertaining to clusters.

"""

import numpy as np
import scipy.cluster, scipy.linalg
import os

from . import fibers, distance, misc, prior
import vtk

def spectralClustering(fiberData, scalarDataList=[], scalarTypeList=[],
                       scalarWeightList=[], k_clusters=50, sigma=[10],
                       n_jobs=-1, dirpath=None, verbose=0):
        """
        Clustering of fibers based on pairwise fiber similarity.
        See paper: "A tutorial on spectral clustering" (von Luxburg, 2007)

        If no scalar data provided, clustering performed based on geometry.
        First element of scalarWeightList should be weight placed for geometry,
        followed by order given in scalarTypeList These weights should sum to
        1.0 (weights given as a decimal value).
        ex.  scalarDataListiberData
             scalarTypeList = [FA, T1]
             scalarWeightList = [Geometry, FA, T1]

        INPUT:
            fiberData - fiber tree of tractography data to be clustered
            scalarDataList - list of scalar data for similarity measurements
            scalarTypeList - list of scalar type for similarity measurements
            scalarWeightList - list of weights for scalar measurements
            k_clusters - number of clusters via k-means clustering
            sigma - width of Gaussian kernel; adjust to alter sensitivity
            n_jobs - number of processes/threads (defaults to use all available
                     resources)
            dirpath - directory to store files
            verbose - verbosity of function

        OUTPUT:
            outputPolydata - polydata containing information from clustering
            clusterIdx - cluster labels of fibers
            fiberData - tree with spatial and quantitative info of fibers
            rejIdx - indices of fibers to reject
        """
        if dirpath is None:
            dirpath = os.getcwd()

        if fiberData.no_of_fibers == 0:
            raise ValueError("Input has 0 fibers!")

        misc.vprint("Starting clustering...", verbose)
        misc.vprint("No. of fibers: %d" % int(fiberData.no_of_fibers), verbose)
        misc.vprint("No. of clusters: %d" % int(k_clusters), verbose)

        # 1. Compute similarty matrix
        W = _pairwiseWeightedSimilarity(fiberData, scalarTypeList,
                                        scalarWeightList, sigma, n_jobs)

        # Outlier detection
        W, rejIdx = _outlierSimDetection(W)

        # 2. Compute degree matrix
        D = _degreeMatrix(W)

        # 3. Compute unnormalized Laplacian
        L = D - W
        del W

        # 4. Compute normalized Laplacian (random-walk)
        D = np.diag(np.divide(1, np.sqrt(np.sum(D, axis=1))))
        Lsym = np.linalg.multi_dot([D, L, D])
        del D, L

        # 5. Compute eigenvalues and eigenvectors of generalized eigenproblem
        # Sort by ascending eigenvalue
        eigval, eigvec = scipy.linalg.eigh(Lsym)
        idx = eigval.argsort()
        eigval, eigvec = eigval[idx], eigvec[:, idx]
        misc.saveEig(dirpath, eigval, eigvec)
        del Lsym, idx

        # 6. Find optimal eigengap and select embedding vector
        gap_idx = _eiggap(eigval)
        emvec = eigvec[:, 1:gap_idx + 1]

        if (k_clusters < gap_idx + 1):
            misc.vprint("WARNING: k-clusters chosen may produce undesirable results",
                        verbose)

        del eigval, eigvec, gap_idx

        # 7. Find clusters using K-means clustering
        centroids, clusterIdx = scipy.cluster.vq.kmeans2(emvec, k_clusters,
                                                         iter=100,
                                                         minit='points')
        centroids, clusterIdx = _sortLabel(centroids, clusterIdx)
        colour = _cluster_to_rgb(centroids)

        misc.vprint("Finished computing clusters...", verbose)

        # 8. Return results
        # Create model with user / default number of chosen samples along fiber
        outputData = fiberData.convertToVTK(rejIdx)
        outputPolydata = _format_outputVTK(outputData, clusterIdx, colour,
                                           centroids)

        # 9. Also add measurements from those used to cluster
        for i in range(len(scalarTypeList)):
            outputPolydata = addScalarToVTK(outputPolydata, fiberData,
                                            scalarTypeList[i], rejIdx=rejIdx)

        return outputPolydata, clusterIdx, fiberData, rejIdx

def spectralPriorCluster(fiberData, priorVTK, templateFlag=False,
                         scalarDataList=[], scalarTypeList=[],
                         scalarWeightList=[], sigma=[10], pflag=True,
                         n_jobs=-1, dirpath=None, verbose=0):
        """
        Clustering of fibers based on pairwise fiber similarity using
        previously clustered fibers via a Nystrom-like method.
        See paper: "A tutorial on spectral clustering" (von Luxburg, 2007)
                   "Spectral grouping using the Nystrom method" (Fowles et al.,
                   2004)

        If no scalar data provided, clustering performed based on geometry.
        First element of scalarWeightList should be weight placed for geometry,
        followed by order given in scalarTypeList. These weights should sum to
        1.0 (weights given as a decimal value).
        ex. scalarDataList
            scalarTypeList = [FA, T1]
            scalarWeightList = [Geometry, FA, T1]

        INPUT:
            fiberData - fiber tree containing tractography data to be clustered
            priorVTK - prior polydata file
            scalarDataList - list with scalar data for similarity measurements;
            scalarTypeList - list with scalar type for similarity type
            scalarWeightList - list with weights for similarity measurements;
            sigma - width of Gaussian kernel; adjust to alter sensitivity
            pflag - flag indicating clustering with priors
            n_jobs - number of processes/threads (defaults to use all available
                     resources)
            dirpath - directory to store files
            verbose - verbosity of function

        OUTPUT:
            outputPolydata - polydata containing information from clustering to
                             be written into VTK
            clusterIdx - cluster labels of fibers
            fiberData - tree with spatial and quantitative info of fibers
        """
        if dirpath is None:
            dirpath = os.getcwd()

        priorData, priorCentroids, priorLabels, subsetIdxes = \
            prior.load(priorVTK, templateFlag)

        if (fiberData.no_of_fibers == 0) or (int(priorData.no_of_fibers) == 0):
            raise ValueError("Input data(s) has 0 fibers!")

        k_clusters = len(priorCentroids)

        misc.vprint("Starting clustering...", verbose)
        misc.vprint("No. of fibers: %d" % int(fiberData.no_of_fibers), verbose)
        misc.vprint("No. of clusters: %d" % int(k_clusters), verbose)

        # 1. Compute similarity matrix
        W, labels = _priorWeightedSimilarity(fiberData, priorData,
                                             scalarTypeList, scalarWeightList,
                                             sigma, pflag, n_jobs)

        misc.vprint("Performing outlier removal...", verbose)
        W, rejIdx = _outlierSimDetection(W, labels=labels,
                                         tflag=templateFlag,
                                         subsetIdxes=subsetIdxes)

        # 2. Identify corresponding cluster indices from similarity
        labels = np.delete(labels, rejIdx)
        clusterIdx = priorLabels[labels]
        fiberData.addClusterInfo(clusterIdx, priorCentroids)
        colour = _cluster_to_rgb(priorCentroids)

        del W

        misc.vprint("Finished clustering...", verbose)
        outputData = fiberData.convertToVTK(rejIdx)
        outputPolydata = _format_outputVTK(outputData, clusterIdx, colour,
                                           priorCentroids)

        misc.vprint("Mapping scalar data if applicable...", verbose)
        # 3. Also add measurements from those used to cluster
        for i in range(len(scalarTypeList)):
            outputPolydata = addScalarToVTK(outputPolydata, fiberData,
                                            scalarTypeList[i], rejIdx=rejIdx)

        return outputPolydata, clusterIdx, fiberData, rejIdx

def addScalarToVTK(polyData, fiberTree, scalarType, fidxes=None, rejIdx=[]):
    """
    Add scalar to all polydata points to be converted to .vtk file.
    This function is different from scalars.addScalar, which only considers
    point used in sampling of fiber.

    INPUT:
        polyData - polydata for scalar measurements to be added to
        fiberTree - the tree containing polydata information
        scalarType - type of quantitative measurement to be aded to polydata
        fidxes - array with fiber indices pertaining to scalar data of
        extracted fibers; default none

    OUTPUT:
        polydata - updated polydata with quantitative information
    """

    data = vtk.vtkFloatArray()
    data.SetName(scalarType.split('/', -1)[-1])

    if fidxes is None:
        fidxes = [i for i in range(fiberTree.no_of_fibers)]

        for i in rejIdx:
            del fidxes[i]

        for fidx in fidxes:
            for pidx in range(0, fiberTree.pts_per_fiber):
                scalarValue = fiberTree.fiberTree[fidx][pidx][scalarType]
                data.InsertNextValue(scalarValue)

    else:
        for i in rejIdx:
            if rejIdx > (len(fidxes) -1):
                continue
            else:
                del fidxes[i]

        for fidx in fidxes:
            for pidx in range(0, fiberTree.pts_per_fiber):
                scalarValue = fiberTree.fiberTree[fidx][pidx][scalarType]
                data.InsertNextValue(float(scalarValue))

    del scalarValue

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
        polyData - extracted cluster in polydata format; no information is
                   retained
    """
    fiberTree = fibers.FiberTree()
    fiberTree.convertFromVTK(inputVTK, pts_per_fiber)

    cluster = fiberTree.getFibers(np.where(clusterIdx == label)[0])
    cluster = fibers.convertFromTuple(cluster)
    polyData = cluster.convertToVTK()

    return polyData

def _pairwiseDistance_matrix(fiberTree, n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Used to compute an NxN distance matrix for all fibers (N) in the input data.

    INPUT:
        fiberTree - tree containing spatial and quantitative information of
                    fibers
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        distances - NxN matrix containing distances between
    """

    distances, _ = distance.fiberDistance(fiberTree.getFibers(
        range(fiberTree.no_of_fibers)), n_jobs=n_jobs)

    if np.diag(distances).all() != 0.0:
        raise ValueError("Diagonals in distance matrix are not equal to 0")

    return distances

def _pairwiseSimilarity_matrix(fiberTree, sigma, n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Computes an NxN similarity matrix for all fibers (N) in the input data.

    INPUT:
        fiberTree - tree containing spatial and quantitative information of
                    fibers
        sigma - width of Gaussian kernel; adjust to alter
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        similarity - NxN matrix containing similarity between fibers based on
                     geometry
    """

    similarity = _pairwiseDistance_matrix(fiberTree, n_jobs=n_jobs)
    similarity = distance.gausKernel_similarity(similarity, sigma)

    # Sanity check
    if np.diag(similarity).all() != 1.0:
        raise ValueError("Diagonals in similarity matrix are not equal to 1")

    return similarity

def _pairwiseQDistance_matrix(fiberTree, scalarType, n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Computes the "pairwise distance" between quantitative points along a fiber.

    INPUT:
        fiberTree - tree containing spatial and quantitative information of fibers
        scalarType - type of quantitative measurements to be used for computation
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        qDistances - NxN matrix containing pairwise distances between fibers
    """

    qDistances = distance.scalarDistance(fiberTree.getScalars(
        range(fiberTree.no_of_fibers), scalarType), n_jobs=n_jobs)

    if np.diag(qDistances).all() != 0.0:
        raise ValueError("Diagonals in distance matrix are not equal to 0")

    return qDistances

def _pairwiseQSimilarity_matrix(fiberTree, scalarType, sigma, n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Computes the similarity between quantitative points along a fiber.

    INPUT:
        fiberTree - tree containing spatial and quantitative information of
                    fibers
        scalarType - type of quantitative measurements to be used for
                     computation
        sigma - width of Gaussian kernel; adjust to alter sensitivity
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        qSimilarity - NxN matrix containing similarity of quantitative
                      measurements between fibers
    """

    qSimilarity = _pairwiseQDistance_matrix(fiberTree, scalarType,
                                            n_jobs=n_jobs)
    qSimilarity = distance.gausKernel_similarity(qSimilarity, sigma)

    # Sanity check
    if np.diag(qSimilarity).all() != 1.0:
        raise ValueError("Diagonals in similarity marix are not equal to 1")

    return qSimilarity

def _priorDistance_matrix(fiberTree, priorTree, pflag=True, n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Used to compute an distance matrix for all fibers (N) in the input data
    through comparison with previously clustered data

    INPUT:
        fiberTree - tree containing spatial and quantitative information of
                    fibers
        priorTree - tree containing spatial and quantitative info from prev.
                    clustered fibers
        pflag - flag to indicate if clustering is performed with priors
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        distances - matrix containing distances between fibers
        labels - corresponding labels of most similar streamlines
    """

    distances, labels = distance.fiberDistance(fiberTree.getFibers(
        range(fiberTree.no_of_fibers)), priorTree.getFibers(
        range(priorTree.no_of_fibers)), pflag=pflag, n_jobs=n_jobs)

    return distances, labels

def _priorSimilarity_matrix(fiberTree, priorTree, sigma, pflag=True, n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Computes a similarity matrix for all fibers (N) in the input data to
    previously clustered fibers

    INPUT:
        fiberTree - tree containing spatial and quantitative information of
                    fibers
        priorTree - tree containing spatial and quantitative info from
                    previously clustered fibers
        sigma - width of Gaussian kernel; adjust to alter
        pflag - flag to indicate if clustering is performed with priors
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        similarity - matrix containing similarity between fibers based on
                     geometry
        labels - corresponding labels of most similar streamlines
    """

    similarities, labels = _priorDistance_matrix(fiberTree, priorTree,
                                                 pflag=pflag, n_jobs=n_jobs)
    similarities = distance.gausKernel_similarity(similarities, sigma)

    return similarities, labels

def _priorQDistance_matrix(fiberTree, priorTree, scalarType, pflag=True,
                           n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Computes the "pairwise distance" between quantitative points along a fiber
    and previously clustered fibers

    INPUT:
        fiberTree - tree containing spatial and quantitative information of
                    fibers
        priorTree - tree containing information on previously clustered fibers
        scalarType - type of quantitative measurements to be used for
                     computation
        pflag - flag to indicate if clustering is performed with priors
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        qDistances - matrix containing pairwise distances between fibers
        labels - corresponding labels of most similar streamlines
    """

    qDistances, labels = distance.scalarDistance(priorTree.getScalars(
        range(priorTree.no_of_fibers), scalarType), pflag=True, n_jobs=n_jobs)
    qDistances = np.array(qDistances)

    return qDistances, labels

def _priorQSimilarity_matrix(fiberTree, priorTree, scalarType, sigma,
                             pflag=True, n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Computes the similarity between quantitative points along a fiber and
    previously clustered fibers

    INPUT:
        fiberTree - tree containing spatial and quantitative information of
                    fibers
        priorsTree - tree containing information on previously clustered fibers
        scalarType - type of quantitative measurements to be used for
                     computation
        sigma - width of Gaussian kernel; adjust to alter sensitivity
        pflag - flag to indicate if clustering is performed with priors
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        qSimilarity - matrix containing similarity of quantitative measurements
        between fibers
        labels - correspondign labels of most similar streamlines
    """

    qDistances, labels = _priorQDistance_matrix(fiberTree, priorTree, scalarType,
                                        pflag=pflag, n_jobs=n_jobs)

    qSimilarity = distance.gausKernel_similarity(qDistances, sigma)

    # Unused variables
    del qDistances

    return qSimilarity, labels

def _degreeMatrix(inputMatrix):
    """ *INTERNAL FUNCTION*
    Computes the degree matrix, D.

    INPUT:
        inputMatrix - adjacency matrix to be used for computation

    OUTPUT:
        degMat - degree matrix to be used to compute Laplacian matrix
    """

    # Determine the degree matrix
    degMat = np.diag(np.sum(inputMatrix, axis=1))

    return degMat

def _cluster_to_rgb(data):
    """ *INTERNAL FUNCTION*
    Generate cluster color from first three components of data

    INPUT:
        data - information used to calculate RGB colours; typically
               eigenvectors or centroids are used

    OUTPUT:
        colour - array containing the RGB values to colour clusters
    """

    colour = data[:, 0:3]

    # Normalize color
    colourMag = np.sqrt(np.sum(np.square(colour), 1))
    colour = np.divide(colour.T, colourMag).T

    # Convert range from 0 to 255
    colour = 127.5 + (colour * 127.5)

    # Unused variable
    del colourMag

    return colour.astype('int')

def _format_outputVTK(polyData, clusterIdx, colour, centroids, rejIdx=[]):
    """ *INTERNAL FUNCTION*
    Formats polydata with cluster index and colour.

    INPUT:
        polyData - polydata for information to be applied to
        clusterIdx - cluster indices to be applied to each fiber within the
                     polydata model
        colour - colours to be applied to each fiber within the polydata model
        centroid - centroid location to associated with each cluster
        rejIdx - indices to reject

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

    for fidx in range(len(clusterIdx)):
        if fidx in rejIdx:
            continue
        dataColour.InsertNextTuple3(colour[clusterIdx[fidx], 0],
                                    colour[clusterIdx[fidx], 1],
                                    colour[clusterIdx[fidx], 2])
        clusterLabel.InsertNextTuple1(int(clusterIdx[fidx]))
        centroid.InsertNextTuple(centroids[clusterIdx[fidx], :])

    polyData.GetCellData().AddArray(dataColour)
    polyData.GetCellData().AddArray(clusterLabel)
    polyData.GetCellData().AddArray(centroid)

    return polyData

def _pairwiseWeightedSimilarity(fiberTree, scalarTypeList=[],
                                scalarWeightList=[], sigma=[10], n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Computes and returns a single weighted similarity matrix.
    Weight list should include weight for distance and sum to 1.

    INPUT:
        fiberTree - tree containing scalar data for similarity measurements
        scalarTypeList - list of scalar type for similarity
                         measurements
        scalarWeightList - list of weights for similarity measurements
        sigma - width of Gaussian kernel; adjust to alter sensitivity
        n_jobs - number of processes/threads (defaults to use all available
                 resources)

    OUTPUT:
        wSimilarity - matrix containing the computed weighted similarity
    """

    if ((scalarWeightList == []) and (scalarTypeList != [])):
        raise ValueError("No weights given for provided measurements!")

    elif ((scalarWeightList != []) and (scalarTypeList == [])):
        raise ValueError("Please also specify measurement(s) type!")

    elif ((((scalarWeightList == [])) and ((scalarTypeList == []))) or
    (scalarWeightList[0] == 1)):
        print("\nCalculating similarity based on geometry.")
        wSimilarity = _pairwiseSimilarity_matrix(fiberTree, sigma[0],
                                                 n_jobs=n_jobs)
        print("\nFinished calculating similarity")

    else:   # Calculate weighted similarity
        if np.sum(scalarWeightList) != 1.0:
            raise ValueError("Weights given do not sum 1!")

        wSimilarity = _pairwiseSimilarity_matrix(fiberTree, sigma[0], n_jobs)
        wSimilarity = wSimilarity * scalarWeightList[0]

        for i in range(len(scalarTypeList)):
            similarity = _pairwiseQSimilarity_matrix(fiberTree,
                scalarTypeList[i], sigma[i+1], n_jobs=n_jobs)

            wSimilarity += similarity * scalarWeightList[i+1]

        del similarity

    if np.diag(wSimilarity).all() != 1.0:
        raise ValueError("Diagonals of weighted similarity are not equal to 1")

    return wSimilarity

def _priorWeightedSimilarity(fiberTree, priorTree, scalarTypeList=[],
                             scalarWeightList=[], sigma=[10], pflag=True,
                             n_jobs=-1):
    """ *INTERNAL FUNCTION*
    Computes and returns a single weighted similarity matrix.
    Weight list should include weight for distance and sum to 1.

    INPUT:
        fiberTree - tree containing scalar data for similarity measurements
        priorTree - tree containing previously clustered tract information
        scalarTypeList - list of scalar type for similarity measurements
        scalarWeightList - list of weights for similarity measurements
        sigma - width of Gaussian kernel; adjust to alter sensitivity
        pflag - flag indicating clustering with priors
        n_jobs - number of processes/threads (defaults to use all available
                 resources)
    OUTPUT:
        wSimilarity - matrix containing the computed weighted similarity
        labels - corresponding labels of most similar streamlines
    """

    if ((scalarWeightList == []) and (scalarTypeList != [])):
        raise ValueError("No weights given for provided measurements!")

    elif ((scalarWeightList != []) and (scalarTypeList == [])):
        raise ValueError("Please also specify measurement(s) type!")

    elif (((scalarWeightList == []) and (scalarTypeList == [])) or
    (scalarWeightList[0] == 1)):
        print("\nCalculating similarity based on geometry.")
        wSimilarity, labels = _priorSimilarity_matrix(fiberTree, priorTree,
                                                      sigma[0], pflag=pflag,
                                                      n_jobs=n_jobs)
        print("\nFinished calculating similarity")

    else:   # Calculate weighted similarity
    # NEEDS TO BE FIXED TO HANDLE HIGHER STREAMLINE COUNTS

        if np.sum(scalarWeightList) != 1.0:
            raise ValueError("Weights given do not sum 1.")

        wSimilarity = _priorSimilarity_matrix(fiberTree, priorTree, sigma[0],
                                              n_jobs=n_jobs)
        wSimilarity = wSimilarity * scalarWeightList[0]

        for i in range(len(scalarTypeList)):
            similarity = _priorQSimilarity_matrix(fiberTree, priorTree,
                scalarTypeList[i], sigma[i+1], n_jobs=n_jobs)

            wSimilarity += similarity * scalarWeightList[i+1]

        del similarity

    return wSimilarity, labels

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

def _outlierSimDetection(W, labels=None, tflag=False, subsetIdxes=None):
    """ * INTERNAL FUNCTION *
    Look for outliers in fibers to reject

    INPUT
        W - similarity matrix
        labels - indices of most similar streamlines for grouping
        tflag - flag for indicating subsetting
        subsetIdxes - subset of indices for subsetting; used with tflag

    OUTPUT:
        W - similarity matrix with removed outliers
        rejIdx - indices of fibers considered outliers
    """

    # Reject fibers that is 1 standard deviations from mean
    if labels is None:
        if tflag is False:
            W_rowsum = np.nansum(W, 1)
            W_outlierthr = np.mean(W_rowsum) - 1 * np.std(W_rowsum)

            rejIdx = np.where(W_rowsum < W_outlierthr)

            # Remove outliers from matrix
            W = np.delete(W, rejIdx[0], axis=0)
            W = np.delete(W, rejIdx[0], axis=1)

            return W, rejIdx
        else:
            rejIdx = [i for i in range(W.shape[0]) if i not in subsetIdxes]
            W = np.delete(W, rejIdx, axis=0)

            return W, rejIdx[0]

    else:
        rejIdx = []
        for label in np.unique(labels):
            temp = np.where(labels == label)
            W_outlierthr = np.mean(W[temp]) - np.std(W[temp])
            for idx in temp[0]:
                if W[idx] < W_outlierthr:
                    rejIdx.append(idx)
                else:
                    continue

        rejIdx = list(np.unique(rejIdx))
        rejIdx.reverse()

        W = np.delete(W, rejIdx)

        return W, rejIdx

def _eiggap(eigval):
    """ * INTERNAL FUNCTION*
    Automatically identify eigengap for stable embedding vector

    Adapted from: https://github.com/mingmingyang/auto_spectral_clustering

    INPUT:
        eigval - eigenvalues from solving eigenproblem

    OUTPUT:
        gap_idx - vectors required for stable embedding vector
    """
    max_gap, gap_idx = 0, 0
    for i in range(1, eigval.size):
        gap = eigval[i] - eigval[i-1]
        if gap > max_gap:
            max_gap = gap
            gap_idx = i - 1
    del gap, max_gap

    return gap_idx
