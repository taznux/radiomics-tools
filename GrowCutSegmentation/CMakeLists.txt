project(FastGrowCut)

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

add_executable(TestSegmentation TestSegmentation.cxx fibheap.cpp)

SET (BUILD_SHARED_LIBS OFF)
SET (CMAKE_VERBOSE_MAKEFILE OFF)

TARGET_LINK_LIBRARIES(TestSegmentation ${ITK_LIBRARIES})
