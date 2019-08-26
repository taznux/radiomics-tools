project(GrowCutSegmentation)

set (EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

include_directories(../ITKUtils)

add_executable(NoduleSegmentation
    NoduleSegmentation.cxx
)

set (BUILD_SHARED_LIBS OFF)
set (CMAKE_VERBOSE_MAKEFILE OFF)

target_link_libraries(NoduleSegmentation ${ITK_LIBRARIES} ITKUtils)


if(codecov_FOUND)
  add_coverage(NoduleSegmentation)
endif()
