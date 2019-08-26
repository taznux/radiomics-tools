project (ITKUtils)

set (LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

add_library(ITKUtils
  ITKUtils.cxx
  itkMaskedSTAPLEImageFilter.hxx
  itkMaskedSpatialObjectToImageFilter.hxx
  itkGrowCutSegmentationImageFilter.hxx
  itkFastGrowCutSegmentationImageFilter.hxx
  fibheap.cpp
)

set (BUILD_SHARED_LIBS OFF)
set (CMAKE_VERBOSE_MAKEFILE OFF)


if(codecov_FOUND)
  add_coverage(ITKUtils)
endif()
