from setuptools import setup
import os

def fullsplit(path, result=None):
    """
    Split a pathname into components (the opposite of os.path.join) in a
    platform-neutral way.
    """
    if result is None:
        result = []
    head, tail = os.path.split(path)
    if head == "":
        return [tail] + result
    if head == path:
        return result
    return fullsplit(head, [tail] + result)

package_dir='qradiomics'

packages = []
for dirpath, dirnames, filenames in os.walk(package_dir):
    # ignore dirnames that start with '.'
    for i, dirname in enumerate(dirnames):
        if dirname.startswith("."):
            del dirnames[i]
    if "__init__.py" in filenames:
        packages.append(".".join(fullsplit(dirpath)))

#print(packages)

setup(name='qradiomics',
        version='0.1',
        packages=packages,
        author='Wookjin Choi',
        install_requires=[
            'setuptools',
            'pandas',
            'scipy',
            'numpy',
            'ipython',
            'matplotlib',
            'ruffus',
            'SimpleITK',
        ]
)
