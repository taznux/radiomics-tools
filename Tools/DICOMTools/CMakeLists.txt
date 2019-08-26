project (DICOMTools)

set (EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

find_package(Boost)
if(Boost_FOUND)
  include_directories( ${Boost_INCLUDE_DIR} )
else()
  message(FATAL_ERROR
    "Cannot build Application without Boost. Please set Boost_DIR.")
endif()

include_directories(../ITKUtils)

add_executable(DICOM2NRRDConverter
    DICOM2NRRDConverter.cxx
)
add_executable(DICOM2NRRDConverterHCC
    DICOM2NRRDConverterHCC.cxx
)
add_executable(DICOM2NRRDConverterNLST
    DICOM2NRRDConverterNLST.cxx
)
add_executable("DICOM-RT2NRRDConverter"
    DICOM-RT2NRRDConverter.cxx
)
add_executable("DICOMTagReader"
    DICOMTagReader.cxx
)


set (BUILD_SHARED_LIBS OFF)
set (CMAKE_VERBOSE_MAKEFILE OFF)

target_link_libraries(DICOM2NRRDConverter ${ITK_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(DICOM2NRRDConverterHCC ${ITK_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(DICOM2NRRDConverterNLST ${ITK_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(DICOM-RT2NRRDConverter ${ITK_LIBRARIES})
target_link_libraries(DICOMTagReader ${ITK_LIBRARIES})



if(codecov_FOUND)
  add_coverage(DICOM2NRRDConverter)
  add_coverage(DICOM2NRRDConverterHCC)
  add_coverage(DICOM2NRRDConverterNLST)
  add_coverage(DICOM-RT2NRRDConverter)
  add_coverage(DICOMTagReader)
endif()
