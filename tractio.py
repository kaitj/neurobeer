""" tractio.py

Module provides input and output functionality for VTK tractography files.

"""

import os
import vtk

VERBOSE = 0

def readVTK(vtkfile):
    """
    Reads tractography .vtk file (vtkPolyData).

    INPUT:
        vtkfile - File of .vtk type containing tractogrpahy polydata

    OUTPUT:
        outData - Polydata stored within the .vtk file
    """

    if VERBOSE:
        print "\nReading", vtkfile, "..."

    filename, ext = os.path.splitext(vtkfile)

    if (ext == '.vtk'):
        vtkReader = vtk.vtkPolyDataReader()
    else:
        print "File format invalid / not recognized."
        return None

    vtkReader.SetFileName(vtkfile)
    vtkReader.Update()
    outData = vtkReader.GetOutput()

    del vtkReader

    if VERBOSE:
        print "Finished reading tractography data..."
        print "Number of fibers found:", outData.GetNumberOfLines()

    return outData

def writeVTK(data, vtkfile):
    """
    Write tractography data to .vtk file (vtkPolyData).

    INPUT:
        data - vtkPolydata to be written to file
        vtkfile - Name of file to be written

    OUTPUT:
        none
    """

    if VERBOSE:
        print "\nWriting", vtkfile, "..."

    filename, ext = os.path.splitext(vtkfile)

    if (ext == '.vtk'):
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileTypeToBinary()
    else:
        print 'Invalid file format'
        return None

    writer.SetFileName(vtkfile)
    writer.SetInputData(data)
    writer.Update()

    del writer

    if VERBOSE:
        print "Finished writing data to ", filename

def readScalar(scalarfile):
    """
    Read input a text file (.txt) containing scalar values pertaining to
    tractography.

    INPUT:
        scalarfile - Text file containing quantitative scalar information

    OUTPUT:
        scalarData - List of scalar values from file
        scalarType - Type of scalar information (ie. FA, T1)
    """

    if VERBOSE:
        print "\nReading", scalarfile, "..."

    scalarType, ext = os.path.splitext(scalarfile)

    if (ext == '.txt'):
        fileReader = open(scalarfile, 'rU')
    else:
        print "File format invalid / not recognized."
        return None

    scalarData = fileReader.readlines()
    fileReader.close()

    for i in range(len(scalarData)):
        scalarData[i] = scalarData[i].rstrip('\n')

    if VERBOSE:
        print "Finished reading scalar data..."

    return scalarData, scalarType
