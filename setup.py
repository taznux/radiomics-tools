from setuptools import setup

setup(name='qradiomcis',
        version='0.1',
        packages=['qradiomics'],
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
