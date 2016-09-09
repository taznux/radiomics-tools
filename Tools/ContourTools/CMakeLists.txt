PROJECT (ContourTools)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

SET (EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

INCLUDE_DIRECTORIES(../ITKUtils)

ADD_EXECUTABLE(ROIGenerator
    ROIGenerator.cxx
)
ADD_EXECUTABLE(ROI2BinImage
    ROI2BinImage.cxx
)
ADD_EXECUTABLE(ROICropImage
    ROICropImage.cxx
)
ADD_EXECUTABLE(ExtractBoundary
    ExtractBoundary.cxx
)
ADD_EXECUTABLE(GTVs2ITV
    GTVs2ITV.cxx
)
ADD_EXECUTABLE(HoleGenerator
    HoleGenerator.cxx
)
ADD_EXECUTABLE(STAPLEComparison
    STAPLEComparison.cxx
)


SET (BUILD_SHARED_LIBS OFF)
SET (CMAKE_VERBOSE_MAKEFILE OFF)

TARGET_LINK_LIBRARIES(ROIGenerator ${ITK_LIBRARIES} ITKUtils)
TARGET_LINK_LIBRARIES(ROICropImage ${ITK_LIBRARIES} ITKUtils)
TARGET_LINK_LIBRARIES(ROI2BinImage ${ITK_LIBRARIES} ITKUtils)
TARGET_LINK_LIBRARIES(ExtractBoundary ${ITK_LIBRARIES} ITKUtils)
TARGET_LINK_LIBRARIES(GTVs2ITV ${ITK_LIBRARIES} ITKUtils)
TARGET_LINK_LIBRARIES(HoleGenerator ${ITK_LIBRARIES} ITKUtils)
TARGET_LINK_LIBRARIES(STAPLEComparison ${ITK_LIBRARIES} ITKUtils)

find_package(codecov)
if(codecov_FOUND)
  add_coverage(ROIGenerator)
  add_coverage(ROICropImage)
  add_coverage(ROI2BinImage)
  add_coverage(ExtractBoundary)
  add_coverage(GTVs2ITV)
  add_coverage(HoleGenerator)
  add_coverage(STAPLEComparison)
endif()
