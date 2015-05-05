PROJECT (DICOMReadAndWrite)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

#------------------------------------------------ 
# ITK
#------------------------------------------------ 

SET(ITK_DIR "/data/research/lu_lab/work/wchoi/software/ITK460-build")

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
	INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
	MESSAGE(FATAL_ERROR
		"Cannot build Application without ITK. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(DICOMReadAndWrite
    DICOMReadAndWrite.cxx
)

SET (BUILD_SHARED_LIBS OFF)
SET (CMAKE_VERBOSE_MAKEFILE OFF)

TARGET_LINK_LIBRARIES(DICOMReadAndWrite ${ITK_LIBRARIES})
