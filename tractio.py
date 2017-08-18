""" tractio.py

Module provides input and output functionality for VTK tractography files.

"""

import os
import vtk

def readVTK(vtkfile):
    """
    Reads tractography .vtk file (vtkPolyData).

    INPUT:
        vtkfile - File of .vtk type containing tractogrpahy polydata

    OUTPUT:
        outData - Polydata stored within the .vtk file
    """
    print("\nReading", vtkfile, "...")

    filename, ext = os.path.splitext(vtkfile)

    if (ext == '.vtk'):
        vtkReader = vtk.vtkPolyDataReader()
    else:
        print("File format invalid / not recognized.")
        return None

    vtkReader.SetFileName(vtkfile)
    vtkReader.Update()
    outData = vtkReader.GetOutput()

    del vtkReader

    print ("Finished reading tractography data...")
    print ("Number of fibers found:", outData.GetNumberOfLines())

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

    print("\nWriting", vtkfile, "...")

    filename, ext = os.path.splitext(vtkfile)

    if (ext == '.vtk'):
        vtkWriter = vtk.vtkPolyDataWriter()
        vtkWriter.SetFileTypeToBinary()
    else:
        print ("Invalid file format")
        return None

    vtkWriter.SetFileName(vtkfile)
    vtkWriter.SetInputData(data)
    vtkWriter.Update()

    del vtkWriter

    print ("Finished writing data to ", filename)

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

    print ("\nReading", scalarfile, "...")

    scalarType, ext = os.path.splitext(scalarfile)

    if (ext == '.txt'):
        fileReader = open(scalarfile, 'rU')
    else:
        print("File format invalid / not recognized.")
        return None

    scalarData = fileReader.readlines()
    fileReader.close()

    for i in range(len(scalarData)):
        scalarData[i] = scalarData[i].rstrip('\n')

    print("Finished reading scalar data...")

    return scalarData, scalarType
