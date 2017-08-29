# Tractography

Development of data driven Tractography tool

## Dependencies

The tool was developed using Python 2.7.12 with the following packages:
```
* Joblib >= 0.10.4.dev0
* Matplotlib >= 1.5.1
* NumPy >= 1.13.1
* Scikit-learn >= 0.19.0
* SciPy >= 0.19.1
* VTK >= 6.3.0
```

All packages, except for VTK, can be installed using the following command:
```
pip install {PACKAGE_NAME}
```

VTK can be downloaded from the official [VTK website](www.vtk.org) and will need to be compiled
via CMake prior to use.
