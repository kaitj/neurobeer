""" fiber.py

Module containing classes and functions related to processing fiber
information.

"""

import numpy as np
import vtk

class Fiber:
    """ Tractography information pertaining to a single fiber """

    def __init__(self):
        self.x = None
        self.y = None
        self.z = None
        self.pts_per_fiber = None

    def getReverseFiber(self):
        """ Determine the reverse order of the fiber with corresponding
        quantitative information

        INPUT:
            none

        OUTPUT:
            none
        """

        fiber = Fiber()

        fiber.x = self.x[::-1]
        fiber.y = self.y[::-1]
        fiber.z = self.z[::-1]

        fiber.pts_per_fiber = self.pts_per_fiber

        return fiber

class FiberArray:
    """ Tractography information pertaining to a group of fibers """

    def __init__(self):
        # Parameters
        self.pts_per_fiber = 2  # Default fibers to end points only

        # Fiber data
        self.fiberArray_x = None
        self.fiberArray_y = None
        self.fiberArray_z = None

        # Quantitative totals
        self.no_of_fibers = 0

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

    def getFiber(self, fiberIdx):
        """ Extract a single fiber from the group with corresponding data.
        Value returned is of class Fiber.

        INPUT:
            fiberIdx - Index of fiber to be extracted

        OUTPUT
            fiber - Single fiber of class Fiber
        """

        fiber = Fiber()

        # Fiber parmaters
        fiber.pts_per_fiber = self.pts_per_fiber

        # Fiber data
        fiber.x = self.fiberArray_x[fiberIdx, :]
        fiber.y = self.fiberArray_y[fiberIdx, :]
        fiber.z = self.fiberArray_z[fiberIdx, :]

        return fiber

    def getFibers(self, fidxes):
        """ Extracts a subset of fibers corresponding to inputted indices.
        Returned fibers are of class fiberArray.

        INPUT:
            fidxes - Indices of subset of fibers to be extracted

        OUTPUT:
            fibers - Subset of fibers of class fiberArray
        """

        fibers = FiberArray()

        # Fiber parameters
        fibers.no_of_fibers = len(fidxes)
        fibers.pts_per_fiber = self.pts_per_fiber

        # Fiber data
        fibers.fiberArray_x = self.fiberArray_x[fidxes]
        fibers.fiberArray_y = self.fiberArray_y[fidxes]
        fibers.fiberArray_z = self.fiberArray_z[fidxes]

        return fibers

    def convertFromVTK(self, inputVTK, pts_per_fiber=None):
        """ Convert input tractography VTK data to array form

        INPUT:
            inputVTK - Tractography polydata
            pts_per_fiber - Number of points to sample along a fiber

        OUTPUT:
            none
        """
        # Points used to discretize along fiber

        if pts_per_fiber is not None:
            self.pts_per_fiber = pts_per_fiber

        # Determine number of lines (assumes all from tractogrpahy)
        self.no_of_fibers = inputVTK.GetNumberOfLines()

        print("\n<fibers.py> Converting polydata to array representation.")
        print ("Fibers:", self.no_of_fibers)
        print("Points along fiber:", self.pts_per_fiber)

        # Initialize fiber storage array: number of fibers, fiber length
        self.fiberArray_x = np.zeros((self.no_of_fibers,
                                                   self.pts_per_fiber))
        self.fiberArray_y = np.zeros((self.no_of_fibers,
                                                    self.pts_per_fiber))
        self.fiberArray_z = np.zeros((self.no_of_fibers,
                                                    self.pts_per_fiber))

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

                self.fiberArray_x[fidx, pidx] = pt[0]
                self.fiberArray_y[fidx, pidx] = pt[1]
                self.fiberArray_z[fidx, pidx] = pt[2]

                pidx += 1

    def convertToVTK(self, scalarArray, scalarType):
        """ Convert fibers in array form to VTK polydata.

        INPUT:
            scalarArray - Variable containing scalar information pertaining
                          to VTK polydata
            scalarType - Type of quantitative scalar (ie. FA, T1) to be
                         included with the poyldata

        OUTPUT:
            outVTK - Tractography polydata in VTK form
        """

        outVTK = vtk.vtkPolyData()
        outPts = vtk.vtkPoints()
        outFibers = vtk.vtkCellArray()
        outScalars = vtk.vtkFloatArray()

        outFibers.InitTraversal()

        # Get fiber information to convert to VTK form
        for fidx in range(0, self.no_of_fibers):
            ptIds = vtk.vtkIdList()

            for pidx in range(0, self.pts_per_fiber):

                idx = outPts.InsertNextPoint(self.fiberArray_x[fidx, pidx],
                                             self.fiberArray_y[fidx, pidx],
                                             self.fiberArray_z[fidx, pidx])

                ptIds.InsertNextId(idx)

                outScalars.InsertNextValue(float(scalarArray.fiberTree_scalar[fidx][pidx][scalarType]))

            outFibers.InsertNextCell(ptIds)

        # Group data into VTK format
        outVTK.SetLines(outFibers)
        outVTK.SetPoints(outPts)
        outVTK.GetPointData().SetScalars(outScalars)

        return outVTK
