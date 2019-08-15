# NeuroBundle Extraction and Evaluation Resource Documentation

The NeuroBundle Extraction and Evaluation Resource, a data-driven tractography clustering and evaluation tool.

[![Docker container available](https://img.shields.io/badge/docker-kaitj/neurobeer-brightgreen.svg?logo=docker&style=flat)](https://hub.docker.com/r/kaitj/neurobeer/tags/)
[![Zenodo](https://zenodo.org/badge/198719661.svg)](https://zenodo.org/badge/latestdoi/198719661)

## Contents
* [Introduction](#introduction)
    * [Disclaimer](#disclaimer)
    * [Support and communication](#support)
* [Installation](#installation)
    * [Containerized package](#container)
* [Workflow](#workflow)
* [References](#references)


## Introduction <a name="introduction"></a>
The NeuroBundle Extraction and Evaluation Resource aims to cluster white matter tractography, leveraging fiber geometry and optionally quantitative MRI information. A spectral clustering algorithm, as described by <sup>1</sup>von Luxburg (2007), is implemented to identify unique tracts from whole-brain tractography. Filtering for short-ranged, U-shaped fibers is performed through adaptation of the filter as described by <sup>2</sup>O'Halloran et al (2017).

### Disclaimer <a name="disclaimer"></a>
No claims are made regarding the correctness of returned output. Results should be interpreted with caution.

### Support and communication <a name="support"></a>
All bugs, concerns, and requests for features can be requested via the github repository found [here](https://github.com/khanlab/neurobeer/issues).

## Installation <a name="installation"></a>
The development of this project was originally written in Python2 and has since been ported to Python3 as of 2019. This package has been tested with both Python2 and Python3.

In order to use the package's library and/or command line interfaces, the latest version of the NeuroBundle Extraction and Evaluation Resource can be downloaded by from the github repository or via git using the following command:

`git clone https://github.com/khanlab/neurobeer`

To install the depedencies needed to use the tool, run the following command:

`pip install -r requirements.txt`

Finally, to install the package, run the following command:

`python setup.py install`

#### Containerized package <a name="container"></a>

NeuroBundle Extraction and Evaluation Resource is also available as a container via both Docker and Singularity.

To use the Docker container, run the following command:

`docker pull kaitj/neurobeer`

To use the Singularity container, users should pull the container from Docker Hub. To do so, run the following command:

`singularity pull docker://kaitj/neurobeer`

_Note: `sudo` may be required to pull or build container._

## Running tractography tool <a name="runmain"></a>
The NeuroBundle Extraction and Evaluation Resource can be used via command line or as a Python library to be used in Python scripts.

## Workflow <a name="workflow"></a>
To use the tool, tractography streamlines should be generated and converted to both .vtk format (for clustering) and .bfloat format (for tracking of scalar quantities along fiber length via Camino<sup>3</sup>).

Scalar data can be generated using the Camino command `tractstats` with the `-tractstat none` in order to generate a text file with corresponding values at each point.

## References <a name="references"></a>
[1] Von Luxburg, U. A tutorial on spectral clustering. Statistics and computing. Stat Comput. 17(4):395-416.

[2] O'Halloran et al. A Method for U-Fiber Quantification from 7T Diffusion-Weighted MRI Data Tested in Subjects with Non-Lesional Focal Epilepsy. Neuroreport. 28(8):457-461.

[3] Cook, PA, et al. Camino: open-source diffusion-mri reconstruction and processing. Proc Intl Soc Magn Reson Med. 2006, 14, 2759.
