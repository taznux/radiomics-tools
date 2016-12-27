! python3 -m pip install --upgrade pip 
! python3 -m pip install conda
! conda config --set always_yes yes --set changeps1 no
! conda install --yes python=$TRAVIS_PYTHON_VERSION pandas scipy numpy ipython matplotlib
! conda install --yes -c bioconda ruffus
! conda install --yes -c SimpleITK SimpleITK 
#conda install pyqt pyside
