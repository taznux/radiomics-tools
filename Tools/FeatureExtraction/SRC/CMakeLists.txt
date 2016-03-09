PROJECT (LungTexture)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

#------------------------------------------------ 
# ITK
#------------------------------------------------ 

SET(ITK_DIR "/data/research/lu_lab/work/wchoi/software/ITK460-build")
IF(WIN32)
  SET(ITK_DIR "C:/Program Files/ITK/lib/cmake/ITK-4.8")
ENDIF(WIN32)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build Application without ITK. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES(../../ITKUtils/SRC)
FIND_LIBRARY(ITKUTILS_LIB ITKUtils ../../ITKUtils/LIB/ REQUIRED)

ADD_EXECUTABLE(FeatureExtraction
    FeatureExtraction.cxx
)


SET (BUILD_SHARED_LIBS OFF)
SET (CMAKE_VERBOSE_MAKEFILE OFF)

TARGET_LINK_LIBRARIES(FeatureExtraction ${ITK_LIBRARIES} ${ITKUTILS_LIB})
