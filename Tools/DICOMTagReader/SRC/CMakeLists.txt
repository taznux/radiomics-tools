PROJECT (DICOMTagReader)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

#------------------------------------------------ 
# ITK
#------------------------------------------------ 
IF(WIN32)
  SET(ITK_DIR "C:/Program Files/ITK/lib/cmake/ITK-4.8")
ENDIF(WIN32)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
	INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
	MESSAGE(FATAL_ERROR
		"Cannot build Application without ITK. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

#------------------------------------------------ 
# Header files
#------------------------------------------------ 
SET (DICOMTagReader_SOURCES_HEADERS
	)

#------------------------------------------------ 
# CPP files
#------------------------------------------------ 
SET (DICOMTagReader_SOURCES_CPP
DICOMTagReader.cxx
)

#--- Group the generated files ---#
SOURCE_GROUP("Template Files"
	REGULAR_EXPRESSION "(.)*\\.txx$"
	
)

ADD_EXECUTABLE(DICOMTagReader
	${DICOMTagReader_SOURCES_HEADERS} 
	${DICOMTagReader_SOURCES_CPP}
	
)

SET (BUILD_SHARED_LIBS OFF)
SET (CMAKE_VERBOSE_MAKEFILE OFF)

TARGET_LINK_LIBRARIES(DICOMTagReader ${ITK_LIBRARIES})
