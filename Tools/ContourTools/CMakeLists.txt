project (ContourTools)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

set (EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

include_directories(../ITKUtils)

add_executable(ROIGenerator
    ROIGenerator.cxx
)
add_executable(ROI2BinImage
    ROI2BinImage.cxx
)
add_executable(ROICropImage
    ROICropImage.cxx
)
add_executable(ExtractBoundary
    ExtractBoundary.cxx
)
add_executable(GTVs2ITV
    GTVs2ITV.cxx
)
add_executable(HoleGenerator
    HoleGenerator.cxx
)
add_executable(STAPLEComparison
    STAPLEComparison.cxx
)


set (BUILD_SHARED_LIBS OFF)
set (CMAKE_VERBOSE_MAKEFILE OFF)

target_link_libraries(ROIGenerator ${ITK_LIBRARIES} ITKUtils)
target_link_libraries(ROICropImage ${ITK_LIBRARIES} ITKUtils)
target_link_libraries(ROI2BinImage ${ITK_LIBRARIES} ITKUtils)
target_link_libraries(ExtractBoundary ${ITK_LIBRARIES} ITKUtils)
target_link_libraries(GTVs2ITV ${ITK_LIBRARIES} ITKUtils)
target_link_libraries(HoleGenerator ${ITK_LIBRARIES} ITKUtils)
target_link_libraries(STAPLEComparison ${ITK_LIBRARIES} ITKUtils)


if(codecov_FOUND)
  add_coverage(ROIGenerator)
  add_coverage(ROICropImage)
  add_coverage(ROI2BinImage)
  add_coverage(ExtractBoundary)
  add_coverage(GTVs2ITV)
  add_coverage(HoleGenerator)
  add_coverage(STAPLEComparison)
endif()
