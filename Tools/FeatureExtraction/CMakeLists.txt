project (FeatureExtraction)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

include_directories(../ITKUtils)

set (EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

add_executable(FeatureExtraction
    FeatureExtraction.cxx
)


set (BUILD_SHARED_LIBS OFF)
set (CMAKE_VERBOSE_MAKEFILE OFF)

target_link_libraries(FeatureExtraction ${ITK_LIBRARIES} ITKUtils)


if(codecov_FOUND)
  add_coverage(FeatureExtraction)
endif()
