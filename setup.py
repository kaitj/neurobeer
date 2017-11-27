from distutils.core import setup

setup(
    # Project information
    name='neurobeer',
    version='0.1.4',
    description='NeuroBundle Extraction and Evaluation Resource',
    packages=['beer/tractography'],
    scripts=['beer/cli/clusterSingle',
                  'beer/cli/clusterPrior',
                  'beer/cli/clusterUFiber',
                  'beer/cli/clusterUFiberPrior'],

    # Metadata
    author='Jason Kai',
    author_email='tkai@uwo.ca',
    url='https://github.com/kaitj/tractography',
)
