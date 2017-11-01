from distutils.core import setup

setup(
    # Project information
    name='tractography',
    version='0.1.3',
    packages=['tractography', 'cli'],
    package_data={'cli': ['clusterSingle', 'clusterPrior', 'clusterUFiber', 'clusterUFiberPrior']},

    # Metadata
    author='Jason Kai',
    author_email='tkai@uwo.ca',
    description='Automated clustering and evaluation tool for white matter tractography',
    url='https://github.com/kaitj/tractography',
)
