""" fiber.py

Module containing classes and functions related to processing fiber
information.

"""

import numpy as np
import vtk
from collections import defaultdict
from . import misc

def tree():
    """
    Creates tree to store quantitative information.

    INPUT:
        none

    OUTPUT:
        none
    """

    return defaultdict(tree)

def convertFromTuple(fiberTuple):
    """
    Converts fiber data in form of type tuple (from extraction) to fiberTree.
    Output is of class FiberTree

    INPUT:
        fiberTuple - tuple containing fiber information to be converted

    OUTPUT:
        fiberTree - fiber information converted to a tree
    """

    fiberTree = FiberTree()

    fiberTree.no_of_fibers = len(fiberTuple[0])
    fiberTree.pts_per_fiber = len(fiberTuple[0][0])

    for fidx in range(fiberTree.no_of_fibers):
        for pidx in range(fiberTree.pts_per_fiber):
                fiberTree.fiberTree[fidx][pidx]['x'] = fiberTuple[0][fidx][pidx]
                fiberTree.fiberTree[fidx][pidx]['y'] = fiberTuple[1][fidx][pidx]
                fiberTree.fiberTree[fidx][pidx]['z'] = fiberTuple[2][fidx][pidx]

    return fiberTree

def calcEndPointSep(fiberData, rejIdx):
    """
    Calculates distance between end points

    INPUT:
        fiberData - fiber tree containing tractography information
        rejIdx - indices of outlier

    OUTPUT:
        DArray -  distance between end points
    """
    endpt = fiberData.pts_per_fiber - 1

    DArray = []
    for fidx in range(fiberData.no_of_fibers):
        if fidx in rejIdx:
            continue
        else:
            x1 = fiberData.fiberTree[fidx][0]['x']
            x2 = fiberData.fiberTree[fidx][endpt]['x']
            y1 = fiberData.fiberTree[fidx][0]['y']
            y2 = fiberData.fiberTree[fidx][endpt]['y']
            z1 = fiberData.fiberTree[fidx][0]['z']
            z2 = fiberData.fiberTree[fidx][endpt]['z']

            DArray.append(np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2))

    return DArray

def calcFiberLength(fiberData, rejIdx=[]):
    """
    Calculates the fiber length via arc length
    NOTE: same function as ufiber module without removing any fibers

    INPUT:
        fiberData - fiber tree containing tractography information
        rejIdx - indices of outlier

    OUTPUT:
        LArray - array containing length of fibers
    """
    no_of_pts = fiberData.pts_per_fiber

    if no_of_pts < 2:
        print("Not enough samples to determine length of fiber")
        raise ValueError

    LArray = []

    for fidx in range(fiberData.no_of_fibers):
        L = 0
        if fidx in rejIdx:
            continue
        else:
            for idx in range(1, no_of_pts):
                x1 = fiberData.fiberTree[fidx][idx]['x']
                x2 = fiberData.fiberTree[fidx][idx - 1]['x']
                y1 = fiberData.fiberTree[fidx][idx]['y']
                y2 = fiberData.fiberTree[fidx][idx - 1]['y']
                z1 = fiberData.fiberTree[fidx][idx]['z']
                z2 = fiberData.fiberTree[fidx][idx - 1]['z']

                L = L + np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

            LArray.append(L)

    return LArray

def addLDRatio(DArray, LArray, polyData):
    """
    Calculates and adds LD Ratio to VTK

    INPUT:
        DArray - array of distances between end points for fibers
        LArray - array of lengths of fibers
        polyData - tractography data to add L/D ratio

    OUTPUT:
        none
    """
    LDScalar = vtk.vtkFloatArray()
    LDScalar.SetNumberOfComponents(1)
    LDScalar.SetName('LDRatio')
    LDRatio = np.divide(DArray, LArray)

    for fidx in range(len(LArray)):
        LDScalar.InsertNextTuple1(LDRatio[fidx])

    polyData.GetCellData().AddArray(LDScalar)

class FiberTree:
    """
    Data pertaining to a group of fibers.
    Value returned is of class FiberTree
    """

    def __init__(self):
        self.fiberTree = tree()

        # Info related to fibers
        self.no_of_fibers = None
        self.pts_per_fiber = None

    def _calc_fiber_indices(self, fiberLength, pts_per_fiber):
        """ *INTERNAL FUNCTION*
        Determine indices to traverse data along a fiber.

        Indices include both end points of the fiber plus evenly spaced points
        along the line. Module determines which indices are wanted based on
        fiber length and desired number of points along the length.

        INPUT:
            fiberLength - number of points along a fiber
            pts_per_fiber - number of desired points along fiber

        OUTPUT:
            idxList - corresponding new indices to traverse along fiber
        """

        # Step length between points along fiber
        stepLength = (fiberLength - 1.0) / (pts_per_fiber - 1.0)

        # Output indices along fiber
        idxList = []
        for idx in range(0, pts_per_fiber):
            idxList.append(idx * stepLength)
        return idxList

    def getFiber(self, fiberIdx):
        """
        Extract a single fiber from the group with corresponding data.
        Value returned is of class Fiber.

        INPUT:
            fiberIdx - index of fiber to be extracted

        OUTPUT
            fiber_x - array of "x" spatial component at each sample
            fiber_y - array of "y" spatial component at each sample
            fiber_z - array of "z" spatial component at each sample
        """

        # Fiber data
        fiber_x = np.zeros(self.pts_per_fiber)
        fiber_y = np.zeros(self.pts_per_fiber)
        fiber_z = np.zeros(self.pts_per_fiber)

        for pidx in range(0, self.pts_per_fiber):
            fiber_x[pidx] = float(self.fiberTree[fiberIdx][pidx]['x'])
            fiber_y[pidx] = float(self.fiberTree[fiberIdx][pidx]['y'])
            fiber_z[pidx] = float(self.fiberTree[fiberIdx][pidx]['z'])

        return fiber_x, fiber_y, fiber_z

    def getFibers(self, fidxes, rejIdx=[]):
        """
        Extracts a subset of fibers corresponding to inputted indices.
        Returned fibers are of class fiberArray.

        INPUT:
            fidxes - Indices of subset of fibers to be extracted

        OUTPUT:
            fiberArray_x - array of "x" spatial component at each sample for
                           fiber bundle
            fiberArray_y - array of "y" spatial component at each sample for
                           fiber bundle
            fiberArray_z - array of "z" spatial component at each sample for
                           fiber bundle
        """

        fiberArray_x = np.zeros((len(fidxes)-len(rejIdx), self.pts_per_fiber))
        fiberArray_y = np.zeros((len(fidxes)-len(rejIdx), self.pts_per_fiber))
        fiberArray_z = np.zeros((len(fidxes)-len(rejIdx), self.pts_per_fiber))

        # Fiber data
        idx = 0
        for fidx in fidxes:
            if fidx in rejIdx:
                continue
            else:
                for pidx in range(0, self.pts_per_fiber):
                    fiberArray_x[idx][pidx] = float(self.fiberTree[fidx][pidx]['x'])
                    fiberArray_y[idx][pidx] = float(self.fiberTree[fidx][pidx]['y'])
                    fiberArray_z[idx][pidx] = float(self.fiberTree[fidx][pidx]['z'])

                idx += 1

        return fiberArray_x, fiberArray_y, fiberArray_z

    def addClusterInfo(self, clusterLabels, centroids):
        """
        Add and save cluster label to fiber tree storing tractography data.

        INPUT:
            clusterLabels - array of cluster labels sorted in fiber index order
            centroids - array of centroids associated with fiber clusters

        OUTPUT:
            none
        """

        uniqueLabels = np.unique(clusterLabels, return_counts=False)

        for label in uniqueLabels:
            for fidx in np.where(clusterLabels == label)[0]:
                self.fiberTree[fidx][str(label)] = label
                self.fiberTree['centroid'][label] = centroids[label]

    def copyScalar(self, fiberData, scalarTypeArray, fidxes=[], rejIdx=[]):
        """ * INTERNAL FUNCTION *
        Copies the scalar information from a different fiber tree. To be used
        for removing outliers while retaining information.

        INPUT:
            fiberData - fiberTree to copy data from
            scalarTypeArray - array of scalar types to copy
            fidxes - array of fiber indices to copy
            rejIdx - array of outlier indices to be exclused; defaults to []

        OUTPUT:
            none
        """

        for Type in scalarTypeArray:
            idx = 0

            # Copy scalars of provided fiber indices
            if fidxes != []:
                for fidx in fidxes:
                    if fidx in rejIdx:
                        continue
                    else:
                        for pidx in range(fiberData.pts_per_fiber):
                            self.fiberTree[idx][pidx][Type] = float(fiberData.fiberTree[fidx][pidx][Type])
                        idx += 1
            # Copy scalars of all fibers
            else:
                for fidx in range(fiberData.no_of_fibers):
                    if fidx in rejIdx:
                        continue
                    else:
                        for pidx in range(fiberData.pts_per_fiber):
                            self.fiberTree[idx][pidx][Type] = float(fiberData.fiberTree[fidx][pidx][Type])
                        idx += 1

    def addScalar(self, inputVTK, scalarData, scalarType, pts_per_fiber=20):
        """
        Add scalar information pertaining to tractography. Values are
        stored with a tree format. This function is dynamic and can add new
        quantitative measurements as needed.

        INPUT:
            inputVTK - tractography polydata to extract corresponding indices
            scalarData - list of scalar values to be stored
            scalarType - type of quantitative scalar (ie. FA, T1)
            pts_per_fiber - number of samples to take along fiber

        OUTPUT:
            none
        """
        inputVTK.GetLines().InitTraversal()
        ptIds = vtk.vtkIdList()

        # Loop over all fibers
        for fidx in range(0, self.no_of_fibers):
            inputVTK.GetLines().GetNextCell(ptIds)
            fiberLength = ptIds.GetNumberOfIds()

            # Loop over pts for ea. fiber
            pidx = 0
            for lineIdx in self._calc_fiber_indices(fiberLength, pts_per_fiber):

                # Find point index
                ptidx = ptIds.GetId(int(round(lineIdx)))
                self.fiberTree[fidx][pidx][scalarType] = float(scalarData[ptidx])

                pidx += 1

    def getScalar(self, fidx, scalarType):
        """
        Extracts scalar information of a specified scalarType pertaining to
        a single fiber.

        INPUT:
            fidx - index corresponding to fiber to extract scalar information
            scalarType - type of quantitative scalar (ie. FA, T1)

        OUTPUT:
            scalarList - list of scalar values indexed by point
        """

        scalarList = np.zeros(self.pts_per_fiber)

        for pidx in range(self.pts_per_fiber):
            scalarList[pidx] = float(self.fiberTree[fidx][pidx][scalarType])

        return scalarList

    def getScalars(self, fidxes, scalarType):
        """
        Extracts scalar information of a specified scalarType pertaining to
        a group of fibers.

        INPUT:
            fidxes - indices corresponding to fibers to extract scalar
                     information from
            scalarType - type of quantitative scalar (ie. FA, T1)

        OUTPUT:
            scalarList - list of scalar values indexed by fiber and point
        """

        scalarList = np.zeros((len(fidxes), self.pts_per_fiber))

        idx = 0
        for fidx in fidxes:
            for pidx in range(0, self.pts_per_fiber):
                scalarList[idx][pidx] = \
                    float(self.fiberTree[fidx][pidx][str(scalarType)])

            idx += 1

        return scalarList

    def convertFromVTK(self, inputVTK, pts_per_fiber=20, verbose=0):
        """
        Convert input tractography VTK data to array form

        INPUT:
            inputVTK - tractography polydata
            pts_per_fiber - number of points to sample along a fiber
            verbose - verbosity of function; 1 to print messages to user.

        OUTPUT:
            none
        """

        self.no_of_fibers = inputVTK.GetNumberOfLines()
        self.pts_per_fiber = pts_per_fiber

        misc.vprint("Converting polydata to array representation.", verbose)
        misc.vprint("Fibers: %d" % int(self.no_of_fibers), verbose)
        misc.vprint("Points sampled along fiber: %d" % int(self.pts_per_fiber),
                     verbose)

        # Loop over all fibers
        inputVTK.GetLines().InitTraversal()
        ptIds = vtk.vtkIdList()
        inputPts = inputVTK.GetPoints()

        for fidx in range(0, self.no_of_fibers):
            inputVTK.GetLines().GetNextCell(ptIds)
            fiberLength = ptIds.GetNumberOfIds()

            # Loop over pts for ea. fiber
            pidx = 0
            for lineIdx in self._calc_fiber_indices(fiberLength,
                                                    self.pts_per_fiber):

                # Perform NN interpolation
                ptidx = ptIds.GetId(int(round(lineIdx)))
                pt = inputPts.GetPoint(ptidx)

                self.fiberTree[fidx][pidx]['x'] = pt[0]
                self.fiberTree[fidx][pidx]['y'] = pt[1]
                self.fiberTree[fidx][pidx]['z'] = pt[2]

                pidx += 1

    def convertToVTK(self, rejIdx=[]):
        """
        Convert fibers in array form to VTK polydata.

        INPUT:
            rejIdx - indices of fibers considered outliers; defaults to []

        OUTPUT:
            outVTK - tractography polydata in VTK form
        """

        outVTK = vtk.vtkPolyData()
        outPts = vtk.vtkPoints()
        outFibers = vtk.vtkCellArray()

        outFibers.InitTraversal()

        # Get fiber information to convert to VTK form
        for fidx in range(0, self.no_of_fibers):
            if fidx in rejIdx:
                continue
            ptIds = vtk.vtkIdList()

            for pidx in range(0, self.pts_per_fiber):
                idx = outPts.InsertNextPoint(self.fiberTree[fidx][pidx]['x'],
                                             self.fiberTree[fidx][pidx]['y'],
                                             self.fiberTree[fidx][pidx]['z'])
                ptIds.InsertNextId(idx)

            outFibers.InsertNextCell(ptIds)

        # Group data into VTK format
        outVTK.SetLines(outFibers)
        outVTK.SetPoints(outPts)

        return outVTK
