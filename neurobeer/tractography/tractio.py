""" tractio.py

Module provides input and output functionality for VTK tractography files.

"""

import os
import vtk

def readVTK(vtkfile, verbose=0):
    """
    Reads tractography .vtk file (vtkPolyData).

    INPUT:
        vtkfile - file of .vtk type containing tractogrpahy polydata
        verbose - verbosity of function; defaults 0

    OUTPUT:
        outData - polydata stored within the .vtk file
    """
    if verbose == 1:
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

    if verbose == 1:
        print "Finished reading tractography data..."
        print "Number of fibers found:", outData.GetNumberOfLines()

    return outData

def writeVTK(data, vtkfile, verbose=0):
    """
    Write tractography data to .vtk file (vtkPolyData).

    INPUT:
        data - vtkPolydata to be written to file
        vtkfile - name of file to be written
        verbose - verbosity of function; defaults 0

    OUTPUT:
        none
    """

    if verbose == 1:
        print"Writing", vtkfile, "..."

    filename, ext = os.path.splitext(vtkfile)

    if (ext == '.vtk'):
        vtkWriter = vtk.vtkPolyDataWriter()
        vtkWriter.SetFileTypeToBinary()
    else:
        print "Invalid file format"
        return None

    vtkWriter.SetFileName(vtkfile)
    vtkWriter.SetInputData(data)
    vtkWriter.Update()

    del vtkWriter

    if verbose == 1:
        print "Finished writing data to ", filename, "\n"

def readScalar(scalarfile, verbose=0):
    """
    Read input a text file (.txt) containing scalar values pertaining to
    tractography.

    INPUT:
        scalarfile - text file containing quantitative scalar information
        verbose - verbosity of function; defaults 0

    OUTPUT:
        scalarData - list of scalar values from file
        scalarType - type of scalar information (ie. FA, T1)
    """
    if verbose == 1:
        print "\nReading", scalarfile, "..."

    scalarType, ext = os.path.splitext(scalarfile)
    scalarType = scalarType.split('_', -1)[-1]

    if (ext == '.txt'):
        fileReader = open(scalarfile, 'rU')
    else:
        print "File format invalid / not recognized."
        return None

    scalarData = fileReader.readlines()
    fileReader.close()

    for i in range(len(scalarData)):
        scalarData[i] = scalarData[i].rstrip('\n')

    if verbose == 1:
        print "Finished reading scalar data..."

    return scalarData, scalarType
