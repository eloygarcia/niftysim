PROJECT(nifty_sim)

cmake_minimum_required(VERSION 2.8)

#-----------------------------------------------------------------------------

if(COMMAND cmake_policy)
   cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#-----------------------------------------------------------------------------
# Default build type
IF (NOT DEFINED CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Default build type, optimisation enabled.")
ENDIF (NOT DEFINED CMAKE_BUILD_TYPE)

#-----------------------------------------------------------------------------
# Platform specific flags
IF (MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_SCL_SECURE_NO_WARNINGS")
ENDIF (MSVC)

#-----------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/source)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/xmlParser)

#-----------------------------------------------------------------------------

OPTION(USE_CUDA "Build GPU solver code (requires CUDA)" OFF)
OPTION(USE_VIZ "Build visualisation classes (requires VTK)" OFF)
OPTION(USE_BOOST "Build with Boost support, required CPU parallel execution" OFF)
OPTION(USE_NAMESPACE_STD "Enables using namespace std directive in tledSimulator header. Deprecated and may be disabled by default in future releases." ON)
OPTION(BUILD_DOXYGEN_DOCS "Create Doxygen API documentation (requires doxygen)" OFF)
OPTION(BUILD_NIFTYSIM_TESTS "Build unit tests" OFF)

IF (USE_CUDA) 
  OPTION(USE_GPU_GP_CONTACT "Use general-purpose CUDA contact modelling code (requires a CUDA device with Fermi architecture or newer, CUDA toolkit version >= 4.2). Please set NIFTYSIM_CUDA_{ARCH,CODE} according to your hardware." OFF)   
ENDIF (USE_CUDA)

#-----------------------------------------------------------------------------

IF (NOT USE_NAMESPACE_STD)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_NOT_USING_NAMESPACE_STD_")
ENDIF (NOT USE_NAMESPACE_STD)

#-----------------------------------------------------------------------------

IF (USE_CUDA AND USE_VIZ)
  # The standard FindCUDA.cmake shipping with all versions up to and including 3.1 is broken when used in conjunction with VTK >= 6.0
  SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
ENDIF (USE_CUDA AND USE_VIZ)

#-----------------------------------------------------------------------------

IF(USE_VIZ)
   ## VTK is needed ##
   FIND_PACKAGE(VTK REQUIRED)
   IF(NOT VTK_FOUND)
      MESSAGE(FATAL_ERROR "VTK package not found")
   ENDIF(NOT VTK_FOUND)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_Visualisation_")

   # INCLUDE(${VTK_USE_FILE})    ## VTK<9.0
   include_directories(${VTK_INCLUDE_DIRS})
ENDIF(USE_VIZ)

#-----------------------------------------------------------------------------

IF (USE_BOOST)  
  IF (USE_CUDA) 
    SET(USE_BOOST OFF)
    MESSAGE(WARNING "CUDA and Boost cannot be used at the same time. Boost is disabled.")
  ELSE ()
    SET(BOOST_LIBS thread date_time)
    FIND_PACKAGE(Boost COMPONENTS ${BOOST_LIBS} REQUIRED)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
    LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_USE_BOOST_")
  ENDIF (USE_CUDA)    
ENDIF (USE_BOOST)

#-----------------------------------------------------------------------------

IF (BUILD_NIFTYSIM_TESTS)
  INCLUDE(${CMAKE_ROOT}/Modules/Dart.cmake)
  ENABLE_TESTING()
  SUBDIRS(testing)
ENDIF (BUILD_NIFTYSIM_TESTS)

#-----------------------------------------------------------------------------

IF(USE_CUDA)
  IF (USE_BOOST) 
    SET(USE_BOOST OFF)
  ENDIF (USE_BOOST) 

  FIND_PACKAGE(CUDA REQUIRED)

  IF (CUDA_VERSION VERSION_GREATER "5.0" OR CUDA_VERSION VERSION_EQUAL "5.0")
    SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-D_CUDA_5PLUS_SDK")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CUDA_5PLUS_SDK")     

    FIND_PATH(CUDA_SDK_COMMON_INCLUDE_DIR
      helper_cuda.h
      PATHS ${CUDA_SDK_SEARCH_PATH}
      PATH_SUFFIXES "common/inc"
      DOC "Location of helper_cuda.h"
      NO_DEFAULT_PATH
      )
    MARK_AS_ADVANCED(CUDA_SDK_COMMON_INCLUDE_DIR)

    INCLUDE_DIRECTORIES(${CUDA_SDK_COMMON_INCLUDE_DIR})
    CUDA_INCLUDE_DIRECTORIES(${CUDA_SDK_COMMON_INCLUDE_DIR})
  ELSE ()
    IF (CUDA_VERSION VERSION_LESS "4.0")
      MESSAGE(WARNING "Your CUDA toolkit version is old, some NiftySim features may not be available or working correctly.")
      SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-D_CUDA_3MINUS_SDK")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CUDA_3MINUS_SDK")
    ENDIF (CUDA_VERSION VERSION_LESS "4.0")

    FIND_PATH(CUDA_CUT_INCLUDE_DIR
      cutil.h
      PATHS ${CUDA_SDK_SEARCH_PATH}
      PATH_SUFFIXES "common/inc"
      DOC "Location of cutil.h"
      NO_DEFAULT_PATH
      )
    MARK_AS_ADVANCED(CUDA_CUT_INCLUDE_DIR)

    INCLUDE_DIRECTORIES(${CUDA_CUT_INCLUDE_DIR})
    CUDA_INCLUDE_DIRECTORIES(${CUDA_CUT_INCLUDE_DIR})
  ENDIF (CUDA_VERSION VERSION_GREATER "5.0" OR CUDA_VERSION VERSION_EQUAL "5.0")

  SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-D_GPU_")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GPU_")
  CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/source)
ENDIF(USE_CUDA)

#-----------------------------------------------------------------------------

IF (USE_GPU_GP_CONTACT)
  IF (CUDA_VERSION VERSION_GREATER "4.2") 
    SET(NIFTYSIM_CUDA_ARCH "compute_80" CACHE STRING "GPU architecture")
    SET(NIFTYSIM_CUDA_CODE "sm_80" CACHE STRING "GPU back-end target")    
    SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-gencode=arch=${NIFTYSIM_CUDA_ARCH},code=${NIFTYSIM_CUDA_CODE};-D_USE_THRUST_;-D_USE_CUDA_ATOMICS;-DGPU_GP_CONTACT")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGPU_GP_CONTACT")
  ELSE (CUDA_VERSION VERSION_GREATER "4.2") 
    MESSAGE(FATAL "Sorry your CUDA toolkit is too old, please upgrade before enabling USE_GPU_GP_CONTACT")
  ENDIF (CUDA_VERSION VERSION_GREATER "4.2") 
ENDIF (USE_GPU_GP_CONTACT)

#-----------------------------------------------------------------------------

IF (BUILD_DOXYGEN_DOCS)
  FIND_PACKAGE(Doxygen)
  IF (DOXYGEN_FOUND)
    IF (USE_CUDA) 
      SET(DOXYGEN_DEFINES "${DOXYGEN_DEFINES} _GPU_")
    ELSE (USE_CUDA)
      IF (USE_BOOST)
	SET(DOXYGEN_DEFINES "${DOXYGEN_DEFINES} _USE_BOOST_")
      ENDIF (USE_BOOST)
    ENDIF (USE_CUDA)

    IF (USE_VIZ) 
      SET(DOXYGEN_DEFINES "${DOXYGEN_DEFINES} _Visualisation_")
    ENDIF (USE_VIZ)

    SET(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
    SET(DOC_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs)
    SET(DOC_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/docs)
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.conf.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf @ONLY)
    ADD_CUSTOM_TARGET(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation..." VERBATIM
      )
  ENDIF (DOXYGEN_FOUND)
ENDIF (BUILD_DOXYGEN_DOCS)
#-----------------------------------------------------------------------------

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CMAKE_USED")

SUBDIRS(xmlParser)
SUBDIRS(source)

#-----------------------------------------------------------------------------

IF (WIN32 AND NOT CYGWIN)
  SET(DEF_INSTALL_CMAKE_DIR CMake)
ELSE ()
  SET(DEF_INSTALL_CMAKE_DIR lib/CMake/NiftySim)
ENDIF (WIN32 AND NOT CYGWIN)
SET(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")
MARK_AS_ADVANCED(INSTALL_CMAKE_DIR)
SET(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
MARK_AS_ADVANCED(INSTALL_LIB_DIR)
SET(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
MARK_AS_ADVANCED(INSTALL_BIN_DIR)
SET(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
MARK_AS_ADVANCED(INSTALL_INCLUDE_DIR)

SET(NiftySim_LIBS tled;xmlParser)
IF (USE_VIZ)
  LIST(APPEND NiftySim_LIBS viz)
ENDIF (USE_VIZ)
MARK_AS_ADVANCED(NiftySim_LIBS)

CONFIGURE_FILE(NiftySimConfig.cmake.in "${PROJECT_BINARY_DIR}/NiftySimConfig.cmake" @ONLY)
CONFIGURE_FILE(UseNiftySim.cmake.in "${PROJECT_BINARY_DIR}/UseNiftySim.cmake" @ONLY)

INSTALL(FILES
  "${PROJECT_BINARY_DIR}/NiftySimConfig.cmake"
  "${PROJECT_BINARY_DIR}/UseNiftySim.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
