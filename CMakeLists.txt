cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
#
# project information
#
project("radiomics-tools" LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11)

#------------------------------------------------
# ITK
#------------------------------------------------

find_package(ITK)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
else()
  message(FATAL_ERROR
    "Cannot build Application without ITK. Please set ITK_DIR.")
endif()

# choose new behaviour for CMP0042
# see http://www.cmake.org/cmake/help/v3.0/policy/CMP0042.html for more details
if (POLICY CMP0042)
	cmake_policy(SET CMP0042 NEW)
endif ()


#
# cmake configuration
#
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/externals/CMake-codecov/cmake" ${CMAKE_MODULE_PATH})

enable_testing()


# enable code coverage
if(EXISTS "${CMAKE_SOURCE_DIR}/externals/CMake-codecov/cmake")
  find_package(codecov)
  set(codecov_FOUND True)
else ()
  set(codecov_FOUND False)
endif ()


#
# recurse into subdirectories
#
add_subdirectory(Tools)

# evaluate coverage
if(codecov_FOUND)
  coverage_evaluate()
endif()
