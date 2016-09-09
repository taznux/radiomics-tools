cmake_minimum_required(VERSION 2.8)

#------------------------------------------------
# ITK
#------------------------------------------------

find_package(ITK)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
else(ITK_FOUND)
  message(FATAL_ERROR
    "Cannot build Application without ITK. Please set ITK_DIR.")
endif(ITK_FOUND)


# choose new behaviour for CMP0042
# see http://www.cmake.org/cmake/help/v3.0/policy/CMP0042.html for more details
if (POLICY CMP0042)
	cmake_policy(SET CMP0042 NEW)
endif (POLICY CMP0042)

#
# project information
#
if (${CMAKE_MAJOR_VERSION} EQUAL 3)
	project("radiomics-tools" LANGUAGES CXX)
else ()
	project("radiomics-tools" CXX)
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
