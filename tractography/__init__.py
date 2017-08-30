"""
The package provides a number of modules to be used in the clustering,
extraction, and evaluation of white matter tractography.

Import of modules for tool

"""

dependencies = ['joblib', 'numpy', 'matplotlib', 'scipy', 'sklearn', 'vtk']
modules = ['cluster', 'distance', 'fibers', 'misc', 'scalars', 'stats', 'tractio']

for dep in dependencies:
    try:
        exec('import %s' % dep)
    except ImportError:
        print '%s is not installed. Functionality of this package will be unavailable..' % dep
        print 'Please installl %s to use this package.' % dep

for module in modules:
    exec('import %s' % module)
