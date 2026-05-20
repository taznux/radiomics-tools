from setuptools import setup, find_packages

setup(name='qradiomics',
        version='0.1',
        packages=find_packages(),
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
