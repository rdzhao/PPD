CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Create project
PROJECT ( CavityDetection )

SET ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules/" )
SET ( TBB_ROOT_DIR "/usr/lib/x86_64-linux-gnu" )

# Find required libraries
FIND_PACKAGE ( OpenVDB REQUIRED )
FIND_PACKAGE ( OpenEXR REQUIRED )  
FIND_PACKAGE ( Boost REQUIRED COMPONENTS system thread )
FIND_PACKAGE ( TBB REQUIRED )
FIND_PACKAGE ( CGAL REQUIRED)

FIND_PACKAGE(OpenMP)
IF (OPENMP_FOUND)
    SET ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
    SET ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" )
ENDIF()

INCLUDE(${CGAL_USE_FILE})

# Project source files
SET ( SRC
  main.cpp
  CavityDetection.h
  CavityDetection.cpp
  nanoflann.hpp
  FlannUtils.h
 )
# Create project 
ADD_EXECUTABLE ( CavityDetection
  ${SRC}
  )

# Project libraries 
TARGET_LINK_LIBRARIES ( CavityDetection
  ${OPENVDB_LIBRARIES}
  ${OPENEXR_LIBRARIES}
  ${ILMBASE_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY_DEBUG}
  ${Boost_THREAD_LIBRARY_DEBUG}
  ${Boost_SYSTEM_LIBRARY_RELEASE}
  ${Boost_THREAD_LIBRARY_RELEASE}
  ${TBB_LIBRARIES}  
  ${CGAL_LIBRARIES}
  ${CGAL_3RD_PARTY_LIBRARIES}
  )

# Project includes
INCLUDE_DIRECTORIES ( ${OPENVDB_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES ( ${OPENEXR_INCLUDE_DIR} )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES ( ${TBB_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES ( .. ) # So they can simulate referencing installed headers
INCLUDE_DIRECTORIES ( . ) # So they can simulate referencing installed headers

ADD_DEFINITIONS ( -pthread -fPIC )


