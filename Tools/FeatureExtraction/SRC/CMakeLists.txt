PROJECT (FeatureExtraction)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

#------------------------------------------------ 
# ITK
#------------------------------------------------ 
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
SET (FeatureExtraction_SOURCES_HEADERS
	)

#------------------------------------------------ 
# CPP files
#------------------------------------------------ 
SET (FeatureExtraction_SOURCES_CPP
FeatureExtraction.cxx
)

#--- Group the generated files ---#
SOURCE_GROUP("Template Files"
	REGULAR_EXPRESSION "(.)*\\.txx$"
	
)

ADD_EXECUTABLE(FeatureExtraction
	${FeatureExtraction_SOURCES_HEADERS} 
	${FeatureExtraction_SOURCES_CPP}
	
)

SET (BUILD_SHARED_LIBS OFF)
SET (CMAKE_VERBOSE_MAKEFILE OFF)

TARGET_LINK_LIBRARIES(FeatureExtraction ${ITK_LIBRARIES})
