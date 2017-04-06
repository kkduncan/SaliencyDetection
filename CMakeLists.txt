cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
##
## SaliencyDetector CMakeLists File
##
## Copyright (c) 2016, KDuncan
## All Rights Reserved.
##

# Adjust some policies
if(COMMAND cmake_policy)
	if(POLICY CMP0010)
		# Bad variable reference syntax
		cmake_policy(SET CMP0010 NEW) 
	endif()	
endif()

# The only configuration type we're interested in at the moment
set(CMAKE_CONFIGURATION_TYPES "Release")
project(SaliencyDetector)

# Turn on Solution Folders (Visual Studio)
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON) 


#=======================================#
# Required Original Headers and Sources #
#=======================================#
set(SAL_INC
	${CMAKE_CURRENT_SOURCE_DIR}/EdgeDetector.h
	${CMAKE_CURRENT_SOURCE_DIR}/SaliencyDetector.h	
)

set(SAL_IMPL	
	${CMAKE_CURRENT_SOURCE_DIR}/EdgeDetector.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SaliencyDetector.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SaliencyMain.cpp
)


#===============#
# Filter Groups #
#===============#
source_group(inc FILES ${SAL_INC})
source_group(src FILES ${SAL_IMPL})


#===========================#
# Third Party Includes etc. #
#===========================#
# Find OpenCV
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
link_directories(${OpenCV_LIBRARY_DIRS})
add_definitions(${OpenCV_DEFINITIONS})


#=======#
# Build #
#=======#
add_executable(SaliencyDetector ${SAL_IMPL} ${SAL_INC})


#======#
# Link #
#======#
target_link_libraries(SaliencyDetector ${OpenCV_LIBS})

