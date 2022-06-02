#!/usr/bin/env bash

export NUM_THREADS=100
export ITK_VERSION=v5.1.2
export PROJECT_DIR="$(pwd)"
export SUPERBUILD_DIR=${PROJECT_DIR}/super-build
export ITK_SOURCE_DIR=${SUPERBUILD_DIR}/itk-${ITK_VERSION}
export ITK_BUILD_DIR=${ITK_SOURCE_DIR}-build
export CONDA_DIR=${PROJECT_DIR}/miniconda/


if ! test -d ${SUPERBUILD_DIR}; then mkdir ${SUPERBUILD_DIR}; fi
cd ${SUPERBUILD_DIR}

# submodule update
git submodule init && git submodule update

# miniconda install
set +e
if [[ "$OSTYPE" == "linux-gnu" ]]; then
  wget "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh" -O miniconda_install.sh
elif [[ "$OSTYPE" == "darwin"* ]]; then
  wget "https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh" -O Miniconda_Install.sh
fi
set -e

bash miniconda_install.sh -b -f -p "${CONDA_DIR}"

# Activate the new environment
PATH="${CONDA_DIR}/bin:${PATH}"

# install python modules
#bash ${PROJECT_DIR}/install_modules.sh
source ${CONDA_DIR}/bin/activate
conda config --set always_yes yes --set changeps1 no
conda install --yes --freeze-installed nomkl 
conda install --yes --freeze-installed pandas scipy numpy ipython matplotlib
conda install --yes --freeze-installed -c conda-forge pydicom pynetdicom 
conda install --yes --freeze-installed -c bioconda ruffus
conda install --yes --freeze-installed -c SimpleITK SimpleITK 
conda clean -afy \
    && find ${CONDA_DIR} -follow -type f -name '*.a' -delete \
    && find ${CONDA_DIR} -follow -type f -name '*.pyc' -delete \
    && find ${CONDA_DIR} -follow -type f -name '*.js.map' -delete 

cd ${PROJECT_DIR}/externals/TCIAExplorer
python3 setup.py install
cd -


if ! test -e ${ITK_SOURCE_DIR}/CMakeLists.txt; then rm -fr $ITK_SOURCE_DIR; fi
if ! test -d ${ITK_SOURCE_DIR}; then git clone --branch ${ITK_VERSION} https://github.com/InsightSoftwareConsortium/ITK.git ${ITK_SOURCE_DIR}; fi
if ! test -d ${ITK_BUILD_DIR}; then mkdir ${ITK_BUILD_DIR}; fi
cd ${ITK_BUILD_DIR}
cmake ${ITK_SOURCE_DIR} -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DModule_ITKReview=ON -DModule_LesionSizingToolkit=ON
make --jobs=${NUM_THREADS} --keep-going
cd -

if ! test -d build; then mkdir build; fi
cd build
cmake ${PROJECT_DIR} -DITK_DIR=${ITK_BUILD_DIR} && make --jobs=${NUM_THREADS}
cd ${PROJECT_DIR}
