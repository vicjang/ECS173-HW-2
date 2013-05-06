#This software is property of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  Do not distribute or copy this software without the consent of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  All rights reserved.  Copyright ©  2007-2013  Owen Carmichael, Chris Schwarz, and the Regents of the University of California.

cmake_minimum_required(VERSION 2.5)

PROJECT(IterativeClosestPoint)
INCLUDE_REGULAR_EXPRESSION("^.*$")

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK) 
IF(VTK_FOUND) 
  INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND) 
  MESSAGE(FATAL_ERROR "VTK not found. Please set VTK_DIR.")
ENDIF(VTK_FOUND)

IF(APPLE)
  FIND_PROGRAM(APPLE_REZ Rez /Developer/Tools)
  ADD_CUSTOM_COMMAND(TARGET gui
    POST_BUILD
    COMMAND ${APPLE_REZ}
    ARGS Carbon.r -o gui)
ENDIF(APPLE)

IF(NOT MSVC)
  ADD_DEFINITIONS( -Wno-deprecated )
ENDIF(NOT MSVC)

ADD_EXECUTABLE(IterativeClosestPoint1 IterativeClosestPoint1.cxx )
TARGET_LINK_LIBRARIES(IterativeClosestPoint1 ITKAlgorithms ITKFEM ITKStatistics ITKIO ITKBasicFilters ITKCommon vtkCommon vtkIO vtkHybrid vtkRendering)

ADD_EXECUTABLE(gui gui.cpp)

TARGET_LINK_LIBRARIES(gui vtkImaging vtkRendering vtkHybrid  vtkGraphics vtkCommon) 
