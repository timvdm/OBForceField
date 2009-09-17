# - Try to find OpenCL
# Once done this will define
#  
#  OPENCL_FOUND        - system has OpenCL
#  OPENCL_INCLUDE_DIR  - the CL include directory
#  OPENCL_LIBRARIES    - Link these to use OpenCL
#   
# If you want to use just CL you can use these values
#  OPENCL_cl_LIBRARY   - Path to OpenCL Library
#  
# On OSX default to using the framework version of opencl
# People will have to change the cache values of 
# OPENCL_cl_LIBRARY to use OpenCL with X11 on OSX

IF (WIN32)
  IF (CYGWIN)

    FIND_PATH(OPENCL_INCLUDE_DIR CL/cl.h )

    FIND_LIBRARY(OPENCL_cl_LIBRARY opencl)

  ELSE (CYGWIN)

    IF(BORLAND)
      SET (OPENCL_cl_LIBRARY import32 CACHE STRING "OpenCL library for win32")
    ELSE(BORLAND)
      SET (OPENCL_cl_LIBRARY opencl CACHE STRING "OpenCL library for win32")
    ENDIF(BORLAND)

  ENDIF (CYGWIN)

ELSE (WIN32)

  IF (APPLE)

    FIND_LIBRARY(OPENCL_cl_LIBRARY OpenCL DOC "OpenCL lib for OSX")
    FIND_PATH(OPENCL_INCLUDE_DIR CL/cl.h DOC "Include for OpenCL on OSX")

  ELSE(APPLE)

    # The first line below is to make sure that the proper headers
    # are used on a Linux machine with the NVidia drivers installed.
    # They replace Mesa with NVidia's own library but normally do not
    # install headers and that causes the linking to
    # fail since the compiler finds the Mesa headers but NVidia's library.
    # Make sure the NVIDIA directory comes BEFORE the others.
    #  - Atanas Georgiev <atanas@cs.columbia.edu>

    FIND_PATH(OPENCL_INCLUDE_DIR CL/cl.h
      ~/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc
      /usr/include
      /usr/local/include
    )

    FIND_PATH(OPENCL_LIBRARY_DIR libOpenCL.so
      ~/NVIDIA_GPU_Computing_SDK/OpenCL/common/lib/Linux64
      ~/NVIDIA_GPU_Computing_SDK/OpenCL/common/lib/Linux34
      /usr/lib
      /usr/local/lib
    )

    FIND_LIBRARY(OPENCL_cl_LIBRARY
      NAMES CL OpenCL 
      PATHS ~/NVIDIA_GPU_Computing_SDK/OpenCL/common/lib/Linux64
            ~/NVIDIA_GPU_Computing_SDK/OpenCL/common/lib/Linux32
            /usr/lib
            /usr/local/lib
            /usr/X11R6/lib
    )

    # On Unix OpenCL most certainly always requires X11.
    # Feel free to tighten up these conditions if you don't 
    # think this is always true.
    # It's not true on OSX.

    IF (OPENCL_cl_LIBRARY)
      IF(NOT X11_FOUND)
        INCLUDE(FindX11)
      ENDIF(NOT X11_FOUND)
      IF (X11_FOUND)
        IF (NOT APPLE)
          SET (OPENCL_LIBRARIES ${X11_LIBRARIES})
        ENDIF (NOT APPLE)
      ENDIF (X11_FOUND)
    ENDIF (OPENCL_cl_LIBRARY)

  ENDIF(APPLE)
ENDIF (WIN32)

SET(OPENCL_FOUND False)
IF(OPENCL_cl_LIBRARY)

    SET( OPENCL_LIBRARIES  ${OPENCL_cl_LIBRARY} ${OPENCL_LIBRARIES})
    SET(OPENCL_FOUND True)

ENDIF(OPENCL_cl_LIBRARY)

# This deprecated setting is for backward compatibility with CMake1.4
SET(OPENCL_INCLUDE_PATH ${OPENCL_INCLUDE_DIR})

MARK_AS_ADVANCED(
  OPENCL_INCLUDE_DIR
  OPENCL_cl_LIBRARY
)

message(STATUS "OPENCL_INCLUDE_DIR: ${OPENCL_INCLUDE_DIR}")
message(STATUS "OPENCL_LIBRARY_DIR: ${OPENCL_LIBRARY_DIR}")
message(STATUS "OPENCL_LIBRARIES: ${OPENCL_LIBRARIES}")
