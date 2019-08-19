from distutils.core import setup
from neurobeer._version import __version__

setup(
    # Project information
    name='neurobeer',
    version=__version__,
    description='NeuroBundle Extraction and Evaluation Resource',
    packages=['neurobeer',
              'neurobeer/tractography'],
    scripts=['neurobeer/cli/clusterSingle',
             'neurobeer/cli/clusterPrior',
             'neurobeer/cli/clusterUFiber',
             'neurobeer/cli/clusterUFiberPrior',
             'neurobeer/cli/tractscalar',
             'neurobeer/cli/vtk2nii',
             'neurobeer/cli/xfmData'],

    # Metadata
    author='Jason Kai',
    author_email='tkai@uwo.ca',
    url='https://github.com/khanlab/neurobeer',
)
