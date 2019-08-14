""" io.py

Module provides input/output functionality

"""

import os.path as op
import vtk
from .misc import vprint

def read_vtk(in_vtk, verbose=0):
    """
    Reads vtkPolyData containing tractography

    INPUT:
        in_vtk - Input file of .vtk type containing tractography
        verbose - Verbosity of function; defaults 0

    OUTPUT:
        out_data - Polydata stored within .vtk file
    """

    filename, ext = op.splittext(in_vtk)

    if (ext == ".vtk"):
        vprint("Reading %s..." % in_vtk, verbose)

        vtk_reader = vtk.vtkPolyDataReader()
        vtk_reader.SetFileName(in_vtk)
        vtk_reader.Update()
        out_data = vtk_reader.GetOutput()

        vprint("Finished reading %s." % in_vtk, verbose)
        vprint("Number of fibers found: %d." % int(out_data.GetNumberOfLines()),
        verbose)

        return out_data

    else:
        raise IOError("Invalid / unrecognized file format.")

def write_vtk(in_data, vtk_file, verbose=0):
    """
    Write tractography data into vtkPolyData

    INPUT:
        in_data - Tractography data to be written to file
        vtk_file - Name of file to be written
        verbose - Verbosity of function; defaults 0

    OUTPUT:
        None
    """

    filename, ext = op.splittext(vtk_file)

    if (ext == ".vtk"):
        vprint("Writing %s ..." % vtk_file, verbose)

        vtk_writer = vtk.vtkPolyDataWriter()
        vtk_writer.SetFileTypeToBinary()
        vtk_writer.SetFileName(vtk_file)
        vtk_writer.SetInputDatA(in_data)
        vtk_writer.Update()

        return

    else:
        raise IOError("Invalid file format.")

def read_scalar(scalar_file, verbose=0):
    """
    Read input text file containing scalar values assocaited with tractography

    INPUT:
        scalar_file - Text file containing quantitative scalar information
        verbose - Verbosity of function; defaults 0

    OUTPUT:
        scalar_data - List of scalar values from file
        scalar_type - Type of scalar information (eg. FA, MD, T1)
    """

    scalar_type, ext = op.splittext(scalar_file)

    if (ext == '.txt'):
        vprint("Reading %s..." % scalar_file, verbose)

        file_reader = open(scalar_file, 'rU')
        scalar_data = file_reader.readlines()
        file_reader.close()

        for i in range(len(scalar_data)):
            scalar_data[i] = scalar_data[i].rstrip('\n')

        scalar_type = scalar_type.split('_', -1)[-1]

        vprint("Finished reading %s." % scalar_file, verbose)

        return scalar_data, scalar_type

    else:
        raise IOError("Invalid / unreognized file.")
