PROJECT (DICOM2NRRDConverter)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

#------------------------------------------------ 
# ITK
#------------------------------------------------ 

IF(WIN32)
  SET(ITK_DIR "C:/Program Files/ITK/lib/cmake/ITK-4.8")
  SET(Boost_INCLUDE_DIR "C:/local/boost_1_59_0")
  SET(Boost_LIBRARY_DIR "C:/local/boost_1_59_0/lib")
ENDIF(WIN32)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build Application without ITK. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(Boost)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
ELSE(Boost_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build Application without Boost. Please set Boost_DIR.")
ENDIF(Boost_FOUND)

ADD_EXECUTABLE(DICOM2NRRDConverter
    DICOM2NRRDConverter.cxx
)
ADD_EXECUTABLE(DICOM2NRRDConverterHCC
    DICOM2NRRDConverterHCC.cxx
)
ADD_EXECUTABLE(DICOM2NRRDConverterNLST
    DICOM2NRRDConverterNLST.cxx
)

SET (BUILD_SHARED_LIBS OFF)
SET (CMAKE_VERBOSE_MAKEFILE OFF)

TARGET_LINK_LIBRARIES(DICOM2NRRDConverter ${ITK_LIBRARIES} ${Boost_LIBRARIES})
TARGET_LINK_LIBRARIES(DICOM2NRRDConverterHCC ${ITK_LIBRARIES} ${Boost_LIBRARIES})
TARGET_LINK_LIBRARIES(DICOM2NRRDConverterNLST ${ITK_LIBRARIES} ${Boost_LIBRARIES})
