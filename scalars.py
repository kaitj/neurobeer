""" scalar.py

Module containing classes and functions related to fiber scalar information.

"""

import numpy as np
import vtk
from collections import defaultdict

def tree():
    """ Creates tree to store quantitative information

    INPUT:
        none

    OUTPUT:
        none
    """

    return defaultdict(tree)

class FiberArrayScalar:
    """ Scalar information pertaining to a group of fibers.
    Value returned is of class FiberScalar. """

    def __init__(self):
        self.fiberTree_scalar = tree()

    def _calc_fiber_indices(self, fiberLength, pts_per_fiber):
        """ *INTERNAL FUNCTION*
        Determine indices to traverse data along a fiber.

        Indices include both end points of the fiber plus evenly spaced points
        along the line. Module determines which indices are wanted based on
        fiber length and desired number of points along the length.

        INPUT:
            fiberLength - Number of points along a fiber
            pts_per_length - Number of desired points along fiber

        OUTPUT:
            idxList - Corresponding new indices to traverse along fiber
        """

        # Step length between points along fiber
        stepLength = (fiberLength - 1.0) / (pts_per_fiber - 1.0)

        # Output indices along fiber
        idxList = []
        for idx in range(0, pts_per_fiber):
            idxList.append(idx * stepLength)

        return idxList

    def addScalar(self, inputVTK, fiberArray, scalarData, scalarType):
        """ Add scalar information pertaining to tractography. Values are
        stored with a tree format. This function is dynamic and can add new
        quantitative measurements as needed.

        INPUT:
            inputVTK - Tractography polydata to extract corresponding indices
            fiberArray - Array of fibers to add scalar to
            scalarData - List of scalar values to be stored
            scalarType - Type of quantitative scalar (ie. FA, T1)

        OUTPUT:
            none
        """

        inputVTK.GetLines().InitTraversal()
        ptIds = vtk.vtkIdList()

        # Loop over all fibers
        for fidx in range(0, fiberArray.no_of_fibers):
            inputVTK.GetLines().GetNextCell(ptIds)
            fiberLength = ptIds.GetNumberOfIds()

            # Loop over pts for ea. fiber
            pidx = 0
            for lineIdx in self._calc_fiber_indices(fiberLength,
                                                    fiberArray.pts_per_fiber):

                # Find point index
                ptidx = ptIds.GetId(int(round(lineIdx)))

                self.fiberTree_scalar[fidx][pidx][scalarType] = \
                    scalarData[ptidx]

                pidx += 1

    def getScalar(self, fiberArray, fidx, scalarType):
        """ Extracts scalar information of a specified scalarType pertaining to
        a single fiber

        INPUT:
            fiber - Lone fiber to get scalar information from
            scalarType - Type of quantitative scalar (ie. FA, T1)

        OUTPUT:
            scalarList - List of scalar values indexed by point
        """

        scalarList = np.zeros(fiberArray.pts_per_fiber)

        for pidx in range(0, fiberArray.pts_per_fiber):

            scalarValue = \
                self.fiberTree_scalar[fidx][pidx][scalarType]
            scalarList[pidx] = float(scalarValue)

        return scalarList

    def getScalars(self, fiberArray, scalarType):
        """ Extracts scalar information of a speficied scalarType pertaining to
        a group of fibers

        INPUT:
            fiberArray - Array of fibers to get scalar from
            scalarType - Type of quantitative scalar (ie. FA, T1)

        OUTPUT:
            scalarList - List of scalar values indexed by fiber and point
        """

        scalarList = np.zeros((fiberArray.no_of_fibers,
                                      fiberArray.pts_per_fiber))

        for fidx in range(0, fiberArray.no_of_fibers):
            for pidx in range(0, fiberArray.pts_per_fiber):

                scalarValue = \
                    self.fiberTree_scalar[fidx][pidx][scalarType]
                scalarList[fidx][pidx] = float(scalarValue)

        return scalarList
