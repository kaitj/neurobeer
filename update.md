# Patch Notes

## Content
* [v1.0.2](#v102)
* [v1.0.1](#v101)
* [v1.0.0](#v100)
* [v0.1.4](#v014)
* [v0.1.3](#v013)
* [v0.1.2](#v012)
* [v0.1.1](#v011)
* [v0.1.0](#v010)

### Version 1.0.2 <a name=v102></a>
* Changed np.dot to np.einsum
* Changed comparison of similarity to lower memory use
* Output directory changed
* Unnecessary variables and functions removed
* Data type changed from float (np.float64) to np.float32
* Eliminated duplications of matrices in calculations
* Removed matrix generation
    * Large memory usage
* Fixed prints for Python3 compatability

### Version 1.0.1 <a name=v101></a>
* Updated prior sampling to subsample from prior/atlas
    * Random subsampling of prior fibers
    * Assignment of cluster index based on similarity (previously Nystrom method)
* Removed interpolation of matrices for visualization

### Version 1.0.0 <a name=v100></a>
* Released ready to use v1.0.0 of tractography tool
* Mapped length:distance ratio for all fibers

### Version 0.1.4 <a name=v014></a>
* Restructured hierarchy; consolidated files under one folder
* Moved CLI commands under scripts in setup.py
* Changed outlier detection using <i>a priori</i> information to make use of distance from centroid
* Output cluster information into CSV
* Added CLI to create template tract by combining subject tractography

### Version 0.1.3 <a name=v013></a>
* Added ufibre module
    * Determines if fibre is u-shaped
    * Approximates fiber length via arc length
    * Calculates distance between end points
* Added outlier detection and rejection

### Version 0.1.2 <a name=v012></a>
* Added in functionality to handle previously clustered / atlas tractography
    * Includes addition of a new prior.py module in library
* Fixed error handling
* Added command line interface for clustering with prior
* Fixed embedding (vector quantization) for priors

### Version 0.1.1 <a name=v011></a>
* Added command line interface functionality
* Set up install script to install library and CLI
* Fixed memory issue caused by figures

### Version 0.1.0 <a name=v010></a>
* Initial deployment
* Library only functionality (Python2)
