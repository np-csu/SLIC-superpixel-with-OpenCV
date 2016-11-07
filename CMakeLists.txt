CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(SLIC-superpixel-with-OpenCV)

SET(OpenCV_DIR "C:\\OpenCV2413\\opencv\\build")

SET(BUILD_SHARED_LIBS on)

FIND_PACKAGE(OpenCV REQUIRED)

MESSAGE(STATUS "OpenCV library status:")
MESSAGE(STATUS "    version: ${OpenCV_VERSION}")
MESSAGE(STATUS "    libraries: ${OpenCV_LIBS}")
MESSAGE(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")

INCLUDE_DIRECTORIES(.)

ADD_LIBRARY(slic
	slic.cpp
	slic.h)
	
TARGET_LINK_LIBRARIES(slic ${OpenCV_LIBS})
	
ADD_EXECUTABLE(demo 
	demo.cpp)
	
TARGET_LINK_LIBRARIES(demo slic)

