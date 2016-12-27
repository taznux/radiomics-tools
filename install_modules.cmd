conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda install --yes python=$TRAVIS_PYTHON_VERSION pandas scipy numpy ipython matplotlib
conda install --yes -c bioconda ruffus
conda install --yes -c SimpleITK SimpleITK 
REM conda install pyqt pyside
