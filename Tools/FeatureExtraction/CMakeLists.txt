PROJECT (FeatureExtraction)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

INCLUDE_DIRECTORIES(../ITKUtils)

SET (EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

ADD_EXECUTABLE(FeatureExtraction
    FeatureExtraction.cxx
)


SET (BUILD_SHARED_LIBS OFF)
SET (CMAKE_VERBOSE_MAKEFILE OFF)

TARGET_LINK_LIBRARIES(FeatureExtraction ${ITK_LIBRARIES} ITKUtils)

find_package(codecov)
add_coverage(FeatureExtraction)
