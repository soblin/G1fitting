############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita` degli Studi di Trento                                   #
#      email: enrico.bertolazzi@ing.unitn.it                               #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.5)

SET( TARGET Clothoid )
SET( CMAKE_EXPORT_COMPILE_COMMANDS ON)
SET( EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
SET( LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
PROJECT( ${TARGET} CXX C )

SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
SET( CMAKE_VERBOSE_MAKEFILE            TRUE )
SET( CMAKE_BUILD_TYPE                  Release )

IF( CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC )
  SET( CMAKE_CXX_FLAGS "-Wall -Wcast-align" )
  IF ( APPLE )
    SET( CMAKE_CXX_FLAGS_RELEASE "-O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-O1 -g -DDEBUG " )
  ELSE()
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -fno-common -O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -fno-common -O1 -g -DDEBUG " )
  ENDIF()
ENDIF()

SET( CMAKE_C_FLAGS         ${CMAKE_CXX_FLAGS} )
SET( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

SET( SOURCES )
FILE( GLOB S ./src/*.cc )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND SOURCES ${RF} )
ENDFOREACH (F ${S})
SET( HEADERS )
FILE( GLOB S ./src/*.hh )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})

INCLUDE_DIRECTORIES( src )

#add_library( ${TARGET} STATIC ${SOURCES} ${HEADERS} )

#ADD_EXECUTABLE( test1 src_tests/test1.cc ${HEADERS} )
#TARGET_LINK_LIBRARIES( test1 ${TARGET} )

#ADD_EXECUTABLE( test2 src_tests/test2.cc ${HEADERS} )
#TARGET_LINK_LIBRARIES( test2 ${TARGET} )

#MESSAGE( STATUS "Using ${SSE_FLAGS} extensions")
#MESSAGE( STATUS "C compiler                  = ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}" )
#MESSAGE( STATUS "C++ compiler                = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}" )
#MESSAGE( STATUS "CMAKE_CXX_FLAGS_RELEASE     = ${CMAKE_CXX_FLAGS_RELEASE}" )
#MESSAGE( STATUS "CMAKE_CXX_FLAGS             = ${CMAKE_CXX_FLAGS}" )
#MESSAGE( STATUS "CMAKE_C_FLAGS_RELEASE       = ${CMAKE_C_FLAGS_RELEASE}" )
#MESSAGE( STATUS "CMAKE_C_FLAGS               = ${CMAKE_C_FLAGS}" )
#MESSAGE( STATUS "CMAKE_SYSTEM_NAME           = ${CMAKE_SYSTEM_NAME}" )
#MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR      = ${CMAKE_SYSTEM_PROCESSOR}" )

# soblin
add_library(
  Clothoid
  SHARED
  src/Clothoid.cc
  src/CubicRootsFlocke.cc
  src/Triangle2D.cc
  )

# for Python
set(LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
add_library(
  clothoid2
  SHARED
  src/Clothoid.cc
  src/CubicRootsFlocke.cc
  src/Triangle2D.cc
  src/python2.cc
  )
set_target_properties(clothoid2 PROPERTIES PREFIX "")
target_link_libraries(clothoid2 boost_python-py27)
target_include_directories(clothoid2 PUBLIC /usr/include/python2.7)

add_library(
  clothoid3
  SHARED
  src/Clothoid.cc
  src/CubicRootsFlocke.cc
  src/Triangle2D.cc
  src/python3.cc
  )
set_target_properties(clothoid3 PROPERTIES PREFIX "")
target_link_libraries(clothoid3 boost_python-py35)
target_include_directories(clothoid3 PUBLIC /usr/include/python3.5)
