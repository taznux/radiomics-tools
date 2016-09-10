PROJECT (ITKUtils)

SET (LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

ADD_LIBRARY(ITKUtils
  ITKUtils.cxx
  itkMaskedSTAPLEImageFilter.hxx
  itkMaskedSpatialObjectToImageFilter.hxx
  itkGrowCutSegmentationImageFilter.hxx
  itkFastGrowCutSegmentationImageFilter.hxx
  fibheap.cpp
)

SET (BUILD_SHARED_LIBS OFF)
SET (CMAKE_VERBOSE_MAKEFILE OFF)


if(codecov_FOUND)
  add_coverage(ITKUtils)
endif()
