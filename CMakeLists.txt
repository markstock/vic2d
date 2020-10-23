#
# vic2d CMake config file
#
CMAKE_MINIMUM_REQUIRED ( VERSION 2.6 )
PROJECT ( vic2d )
ENABLE_LANGUAGE (Fortran)

# User-serviceable variables here
SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE "Release")
ENDIF ()
SET (CMAKE_BUILD_TYPE "Release" CACHE STRING "Select which configuration to build")
SET (CMAKE_INSTALL_PREFIX CACHE PATH "Installation location for binaries, sample inputs, and licenses")
SET (USE_OMP FALSE CACHE BOOL "Use OpenMP multithreading")
SET (CMAKE_VERBOSE_MAKEFILE on)
SET (CMAKE_EXPORT_COMPILE_COMMANDS on)

IF (CMAKE_COMPILER_IS_GNUCXX)
  SET (CMAKE_C_FLAGS "-Wall -Wformat -std=c99")
  SET (CMAKE_C_FLAGS_DEBUG "-O0 -funroll-loops -ffast-math -g -ggdb -fbounds-check")
  SET (CMAKE_C_FLAGS_RELEASE "-O2 -funroll-loops -ffast-math -fomit-frame-pointer")
  SET (CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -funroll-loops -ffast-math -fomit-frame-pointer -g -ggdb3")
  SET (CMAKE_Fortran_FLAGS "-Wall -ffree-line-length-none -Wextra -fno-frontend-optimize")
  SET (CMAKE_Fortran_FLAGS_DEBUG "-O0 -funroll-loops -ffast-math -g -ggdb -fbounds-check")
  SET (CMAKE_Fortran_FLAGS_RELEASE "-O2 -funroll-loops -ffast-math -fomit-frame-pointer")
  SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -funroll-loops -ffast-math -fomit-frame-pointer -g -ggdb3")
ELSEIF (CMAKE_C_COMPILER_ID STREQUAL "Clang")
  SET (CMAKE_C_FLAGS "-Wall -Wformat -std=c99")
  SET (CMAKE_C_FLAGS_DEBUG "-O0 -g -ggdb")
  SET (CMAKE_C_FLAGS_RELEASE "-O3")
  SET (CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g -ggdb")
ELSEIF (MSVC)
  SET (CMAKE_C_FLAGS "/Ox")
ENDIF ()

# suppress the FindOpenMP warnings
cmake_policy(SET CMP0054 NEW)

IF ( WIN32 )

  SET ( PLATFORM_LIBS ${FRAMEWORK_LIBS} )

  # OpenMP for multithreading
  IF( USE_OMP )
    FIND_PACKAGE (OpenMP REQUIRED)
    IF (OPENMP_FOUND)
      #MESSAGE ("OPENMP FOUND")
      SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ENDIF ()
  ENDIF()

ELSEIF ( APPLE )

  INCLUDE_DIRECTORIES ( /System/Library/Frameworks )
  FIND_LIBRARY (COCOA_LIBRARY Cocoa)
  FIND_LIBRARY (IOKIT_LIBRARY IOKit)
  MARK_AS_ADVANCED (COCOA_LIBRARY IOKIT_LIBRARY)

  SET (PLATFORM_LIBS ${FRAMEWORK_LIBS} ${COCOA_LIBRARY} ${IOKIT_LIBRARY})

  # OpenMP for multithreading
  IF( USE_OMP )
    FIND_PACKAGE (OpenMP REQUIRED)
    IF (OPENMP_FOUND)
      #MESSAGE ("OPENMP FOUND")
      SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ENDIF ()
  ENDIF()

ELSEIF ( UNIX )

  # OpenMP for multithreading
  IF( USE_OMP )
    FIND_PACKAGE( OpenMP REQUIRED )
    IF( OPENMP_FOUND )
      #message("OPENMP FOUND")
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    ENDIF()
  ELSE()
    # suppress the warning
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas")
    # need to find pthreads (for std::async?)
    FIND_PACKAGE( Threads REQUIRED )
  ENDIF()

  SET( EXTERNAL_LIBS png z )

ENDIF ()

INCLUDE_DIRECTORIES ( "src" )
ADD_DEFINITIONS (${CPREPROCDEFS})

FILE(GLOB HEADERS "src/*.h" )
SET(LIBSRCS2D ${HEADERS}
            "src/libvicmoc2d.c"
            "src/gr2.f"
            "src/utility.c"
            "src/maskops.c"
            "src/mud2sp_extern.c"
            "src/mud2sp_full.f" )
SET(SRCS2D ${HEADERS}
            "src/vic2d.c"
            "src/particles.c"
            "src/inout.c" )
SET(LIBSRCS3D ${HEADERS}
            "src/libvicmoc3d.c"
            "src/gr3.f"
            "src/utility.c"
            "src/mud3sp.f" )

# targets
ADD_LIBRARY ( "vicmoc2d" STATIC ${LIBSRCS2D} )
ADD_EXECUTABLE ( "vic2d" ${SRCS2D} )
TARGET_LINK_LIBRARIES( "vic2d" "vicmoc2d" ${EXTERNAL_LIBS} )

ADD_LIBRARY ( "vicmoc3d" STATIC ${LIBSRCS3D} )
