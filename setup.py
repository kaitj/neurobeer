from distutils.core import setup

setup(
    # Project information
    name='neurobeer',
    version='0.1.4',
    description='NeuroBundle Extraction and Evaluation Resource',
    packages=['neurobeer/tractography'],
    scripts=['neurobeer/cli/clusterSingle',
                  'neurobeer/cli/clusterPrior',
                  'neurobeer/cli/clusterUFiber',
                  'neurobeer/cli/clusterUFiberPrior',
                  'neurobeer/cli/createTemplateTract'],

    # Metadata
    author='Jason Kai',
    author_email='tkai@uwo.ca',
    url='https://github.com/kaitj/neurobeer',
)
