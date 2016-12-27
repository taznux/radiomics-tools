#!/usr/bin/env bash

export NUM_THREADS=4
export ITK_VERSION=v4.10.0
export PROJECT_DIR="$(pwd)"
export SUPERBUILD_DIR=${PROJECT_DIR}/super-build
export ITK_SOURCE_DIR=${SUPERBUILD_DIR}/itk-${ITK_VERSION}
export ITK_BUILD_DIR=${ITK_SOURCE_DIR}-build


if ! test -d ${SUPERBUILD_DIR}; then mkdir ${SUPERBUILD_DIR}; fi
cd ${SUPERBUILD_DIR}

# submodule update
git submodule init && git submodule update

# miniconda install
set +e
if [[ "$OSTYPE" == "linux-gnu" ]]; then
  curl "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh" -o miniconda_install.sh
  if [ $? -ne 0 ]; then
      curl "http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh" -o miniconda_install.sh
  fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
  curl "https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh" -o miniconda_install.sh
  if [ $? -ne 0 ]; then
      curl "http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh" -o miniconda_install.sh
  fi
fi
set -e

bash miniconda_install.sh -b -f -p "${SUPERBUILD_DIR}/miniconda/"

# Activate the new environment
PATH=${SUPERBUILD_DIR}"/miniconda/bin":$PATH

# install python modules
bash ${PROJECT_DIR}/install_modules.sh

conda install cryptography

cd ${PROJECT_DIR}/externals/TCIAExplorer
python3 setup.py install
cd -


if ! test -e ${ITK_SOURCE_DIR}/CMakeLists.txt; then rm -fr $ITK_SOURCE_DIR; fi
if ! test -d ${ITK_SOURCE_DIR}; then git clone --branch ${ITK_VERSION} git://itk.org/ITK.git ${ITK_SOURCE_DIR}; fi
if ! test -d ${ITK_BUILD_DIR}; then mkdir ${ITK_BUILD_DIR}; fi
cd ${ITK_BUILD_DIR}
cmake ${ITK_SOURCE_DIR} -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DModule_ITKReview=ON
make --jobs=${NUM_THREADS} --keep-going
cd -

if ! test -d build; then mkdir build; fi
cd build
cmake ${PROJECT_DIR} -DITK_DIR=${ITK_BUILD_DIR} && make
cd ${PROJECT_DIR}
