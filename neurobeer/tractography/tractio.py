""" tractio.py

Module provides input/output functionality

"""

import os.path as op
import vtk
from . import misc

def readVTK(in_vtk, verbose=0):
    """
    Reads vtkPolyData containing tractography

    INPUT:
        in_vtk - input file of .vtk type containing tractography
        verbose - verbosity of function; defaults 0

    OUTPUT:
        out_data - polydata stored within .vtk file
    """

    filename, ext = op.splitext(in_vtk)

    if (ext == ".vtk"):
        misc.vprint("Reading %s..." % in_vtk, verbose)

        vtk_reader = vtk.vtkPolyDataReader()
        vtk_reader.SetFileName(in_vtk)
        vtk_reader.Update()
        out_data = vtk_reader.GetOutput()

        misc.vprint("Finished reading %s." % in_vtk, verbose)
        misc.vprint("Number of fibers found: %d." % int(out_data.GetNumberOfLines()),
        verbose)

        return out_data

    else:
        raise IOError("Invalid / unrecognized file format.")

def writeVTK(in_data, vtk_file, verbose=0):
    """
    Write tractography data into vtkPolyData

    INPUT:
        in_data - tractography data to be written to file
        vtk_file - name of file to be written
        verbose - verbosity of function; defaults 0

    OUTPUT:
        none
    """

    filename, ext = op.splitext(vtk_file)

    if (ext == ".vtk"):
        misc.vprint("Writing %s ..." % vtk_file, verbose)

        vtk_writer = vtk.vtkPolyDataWriter()
        vtk_writer.SetFileTypeToBinary()
        vtk_writer.SetFileName(vtk_file)
        vtk_writer.SetInputData(in_data)
        vtk_writer.Update()

        return

    else:
        raise IOError("Invalid file format.")

def readScalar(scalar_file, verbose=0):
    """
    Read input text file containing scalar values assocaited with tractography

    INPUT:
        scalar_file - text file containing quantitative scalar information
        verbose - verbosity of function; defaults 0

    OUTPUT:
        scalar_data - list of scalar values from file
        scalar_type - type of scalar information (eg. FA, MD, T1)
    """

    scalar_type, ext = op.splitext(scalar_file)

    if (ext == '.txt'):
        misc.vprint("Reading %s..." % scalar_file, verbose)

        file_reader = open(scalar_file, 'rU')
        scalar_data = file_reader.readlines()
        file_reader.close()

        for i in range(len(scalar_data)):
            scalar_data[i] = scalar_data[i].rstrip('\n')

        scalar_type = scalar_type.split('_', -1)[-1]

        misc.vprint("Finished reading %s." % scalar_file, verbose)

        return scalar_data, scalar_type

    else:
        raise IOError("Invalid / unreognized file.")
