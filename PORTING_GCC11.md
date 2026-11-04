## C++ compatibility: From GCC11 to GCC17

## Current validated environment (2026-11-04)

- OS: Ubuntu 24.04
- C++ compiler: GCC/G++ 13.3.0
- CMake: 3.28.3
- CUDA Toolkit: 12.6.1

Notes:
- This document keeps the original migration history (GCC11/VTK9) and extends it
  with the CUDA modernization that was applied afterwards.
- The current [Dockerfile](./Dockerfile) still reflects an older base
  (Ubuntu 20.04, CUDA 11.7.1, VTK 8.2) and should not be considered the current
  validated runtime baseline.

Patch 001
----------
File: ./source/tledXMLExporter.tpp
Line: 24
Reason: Old stream cast trick not accepted by GCC11.
Fix: Replace temporary ostringstream cast with explicit stream object.
Old version:
    XMLNode tledXMLExporter<TInput>::AddNumericAttribute(XMLNode &r_node, const std::string &name, const TValue value) {
        r_node.addAttribute(name.c_str(), static_cast<std::ostringstream&>(std::ostringstream() << value).str().c_str());
        return r_node;
    }
New version:
    XMLNode tledXMLExporter<TInput>::AddNumericAttribute(XMLNode &r_node, const std::string &name, const TValue value) {
        std::ostringstream oss;
        oss << value;
        r_node.addAttribute(name.c_str(), oss.str().c_str());
        return r_node;
    }

Patch 002
---------
File:./source/tledModel.cpp
Line: 1019
Issue: Potential buffer overflow warning emitted by GCC11 for sprintf("%15.14f").
Fix: Replace sprintf with snprintf and increase buffer size.
Category: GCC11 hardening / safety checks.
Old version:
    void tledModel::SetTimeStep(double dt)
    {
        char* str = new char[15];
        sprintf(str,"%15.14f",dt);
        xModel.getChildNode("SystemParams").getChildNode("TimeStep").updateText(str);
        delete str;
    }

New Version:
    void tledModel::SetTimeStep(double dt)
    {
        char str[64];
        snprintf(str, sizeof(str), "%15.14f", dt);
        xModel.getChildNode("SystemParams").getChildNode("TimeStep").updateText(str);
    }

## VTK compatibility: From VTK8.2 to VTK9.1

Patch 003
---------
File: ./CMakeLists.txt
From Line: 85  to line: 187 (after: SUBDIRS(source))
Issue: NiftySim tests are compiled before the main application.
Fix: Relocate `IF (BUILD_NIFTYSIM_TESTS) .. ENDIF()`.

Patch 004
---------
File: ./CMakeLists.txt
From Line: 110 to line: 128
Reason: Remove compatibility with CUDA<5.0 in `IF(USE_CUDA) .. ENDIF()`.

Patch 005
---------
File: ./CMakeLists.txt
From Line: 138 to line: 147
Reason: Remove error message for CUDA<4.2 in `IF(USE_GPU_GP_CONTACT) .. ENDIF()`.

Patch 006
---------
File: ./CMakeLists.txt
From Line:  56 to line: 66
Reason: Updating CMake files to make them compatible with VTK9
TODO: Update `$VTK_DIR`.

Old code: vtk==8.2
    IF(USE_VIZ)
        FIND_PACKAGE(VTK REQUIRED)
        IF(NOT VTK_FOUND)
            MESSAGE(FATAL_ERROR "VTK package not found")
        ENDIF(NOT VTK_FOUND)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_Visualisation_")
        include_directories(${VTK_INCLUDE_DIRS})
    ENDIF(USE_VIZ)

New code: vtk==9.1
    IF(USE_VIZ)
        find_package(VTK REQUIRED COMPONENTS
            CommonCore
            CommonDataModel
            FiltersCore
            FiltersSources
            InteractionStyle
            RenderingCore
            RenderingOpenGL2
            InteractionStyle
            IOLegacy
            IOImage
        )
        
        ADD_DEFINITIONS(-D_Visualisation_)
        MESSAGE(STATUS "VTK_VERSION=${VTK_VERSION}")

        IF(NOT VTK_FOUND)
            MESSAGE(FATAL_ERROR "VTK package not found")
        ENDIF(NOT VTK_FOUND)
        
        IF(VTK_VERSION VERSION_LESS "9.0")
            INCLUDE(${VTK_USE_FILE})    ## VTK<9.0
            include_directories(${VTK_INCLUDE_DIRS})
        ENDIF()
    ENDIF(USE_VIZ)


Patch 007
---------
File: ./source/CMakeLists.txt
From Line: 427 to line: 474
Reason: Updating CMake files to make them compatible with VTK9
TODO: Update `$VTK_DIR`.

Old code: 
    IF(USE_VIZ)
        SET(NAME viz)
        SET( CPP_SOURCE_FILES_WITH_HEADERS
            tledModelViewer
            tledKeyPressInteractorStyle
            tledVTKCylSource
            tledVTKUltrasoundProbeSource
            tledVTKPlateSource
            tledVTKSurfaceSource
            tledContactSurfaceExporter
            tledVTKMeshSource
        )
        FOREACH( SOURCE_FILE ${CPP_SOURCE_FILES_WITH_HEADERS} )
            SET( VIZ_LIST ${VIZ_LIST} ${SOURCE_FILE}.cpp )
            INSTALL( FILES ${SOURCE_FILE}.h DESTINATION include )
        ENDFOREACH( SOURCE_FILE )

        ADD_LIBRARY( ${NAME} ${VIZ_LIST} )

        IF (VTK_MAJOR_VERSION LESS 6) 
            TARGET_LINK_LIBRARIES(${NAME} vtkHybrid)
        ELSE ()
            TARGET_LINK_LIBRARIES(${NAME} ${VTK_LIBRARIES})
        ENDIF (VTK_MAJOR_VERSION LESS 6) 

        INSTALL(TARGETS ${NAME}
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        )
    ENDIF(USE_VIZ)

Current version:
    IF(USE_VIZ)
        SET(NAME viz)
        SET( CPP_SOURCE_FILES_WITH_HEADERS
            tledModelViewer
            tledKeyPressInteractorStyle
            tledVTKCylSource
            tledVTKUltrasoundProbeSource
            tledVTKPlateSource
            tledVTKSurfaceSource
            tledContactSurfaceExporter
            tledVTKMeshSource
        )
        FOREACH( SOURCE_FILE ${CPP_SOURCE_FILES_WITH_HEADERS} )
            SET( VIZ_LIST ${VIZ_LIST} ${SOURCE_FILE}.cpp )
            INSTALL( FILES ${SOURCE_FILE}.h DESTINATION include )
        ENDFOREACH( SOURCE_FILE )

        ADD_LIBRARY( ${NAME} ${VIZ_LIST} )

        TARGET_LINK_LIBRARIES(${NAME}
            PUBLIC
            ${VTK_LIBRARIES}
        )
        
        target_include_directories(tled
            PUBLIC
            /usr/include/vtk-9.1
            # ${VTK_DIR}
        )

        if(VTK_VERSION VERSION_GREATER_EQUAL 9.0)
            vtk_module_autoinit(
                TARGETS ${NAME}
                MODULES ${VTK_LIBRARIES}
            )
            endif()

        INSTALL(TARGETS ${NAME}
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        )
    ENDIF(USE_VIZ)



Patch 008
---------
File: ./testing/CMakeLists.txt
From Line: 55 to line: 69
Reason: Updating CMake files to make them compatible with VTK9
TODO: Update `$VTK_DIR`.

Old code:
    FOREACH (TEST ${NO_VTK_TESTS}) 
        IF (USE_CUDA)
            CUDA_ADD_EXECUTABLE(${TEST} ${TEST}.cpp)
        ELSE (USE_CUDA)
            ADD_EXECUTABLE(${TEST} ${TEST}.cpp)
        ENDIF (USE_CUDA)
            
        IF (USE_VIZ)
            TARGET_LINK_LIBRARIES(${TEST} unitTests tled viz)
        ELSE (USE_VIZ)
            TARGET_LINK_LIBRARIES(${TEST} unitTests tled)
        ENDIF (USE_VIZ)  
        STRING(REGEX REPLACE "^test" "" TESTNAME ${TEST})
        ADD_TEST(NAME ${TESTNAME} COMMAND ${TEST})
    ENDFOREACH (TEST) 

Current code:
    FOREACH (TEST ${NO_VTK_TESTS}) 
        IF (USE_CUDA)
            CUDA_ADD_EXECUTABLE(${TEST} ${TEST}.cpp)
        ELSE (USE_CUDA)
            ADD_EXECUTABLE(${TEST} ${TEST}.cpp)
        ENDIF (USE_CUDA)
            
        IF (USE_VIZ)
            TARGET_LINK_LIBRARIES(${TEST} unitTests tled viz)
            target_include_directories(unitTests
                PRIVATE
                /usr/include/vtk-9.1
                #${VTK_DIR}
            )
        ELSE (USE_VIZ)
            TARGET_LINK_LIBRARIES(${TEST} unitTests tled)
        ENDIF (USE_VIZ)  
        STRING(REGEX REPLACE "^test" "" TESTNAME ${TEST})
        ADD_TEST(NAME ${TESTNAME} COMMAND ${TEST})
    ENDFOREACH (TEST) 


Patch 009
---------
File: ./testing/testVTKMeshExporter.cpp
Line: 129 & 141
Error:  error: 'class vtkFloatArray' has no member named 'GetTupleValue'
Fix: Replace `GetTupleValue` with `GetTuple`.

Old code:
    p_array->GetTupleValue(vInd, tmpTuple);

New code:
    p_array->GetTypedTuple(vInd, tmpTuple);

Patch 010
---------
File: ./CMakeLists.txt
Line: 66 (Within IF(USE_VIZ) .. ENDIF())
Reason: Update `$VTK_DIR` in CMake files to make them compatible with VTK9.
Include:

 get_target_property(
    VTK_INCLUDE_DIR
    VTK::CommonCore
    INTERFACE_INCLUDE_DIRECTORIES
)

  message(STATUS "VTK_INCLUDE_DIR=${VTK_INCLUDE_DIR}")
  
  set(NIFTYSIM_VTK_INCLUDE_DIR
    ${VTK_INCLUDE_DIR}
    CACHE INTERNAL "VTK include directory"
)

Patch 011
---------
File: ./testing/CMakeLists.txt
Line: 67
Reason: Update `$VTK_DIR` in CMake files to make them compatible with VTK9.
Modify:

    target_include_directories(unitTests
         PUBLIC
         /usr/include/vtk-9.1
    )

by:

    target_include_directories(unitTests
         PUBLIC
         ${NIFTYSIM_VTK_INCLUDE_DIR}
    )

Patch 012
---------
File: ./source/CMakeLists.txt
Line: 454
Reason: Update `$VTK_DIR` in CMake files to make them compatible with VTK9.
Modify:

    target_include_directories(tled
         PUBLIC
         /usr/include/vtk-9.1
    )

by:

    target_include_directories(tled
         PUBLIC
         ${NIFTYSIM_VTK_INCLUDE_DIR}
    )

## CUDA

Patch 013
---------
File: ./CMakeLists.txt
Line: 47
Reason: Remove warning incompatibilities between CUDA and VTK.
Old code:
    IF (USE_CUDA AND USE_VIZ)
        # The standard FindCUDA.cmake shipping with all versions up to and including 3.1 is broken when used in conjunction with VTK >= 6.0
        SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
    ENDIF (USE_CUDA AND USE_VIZ)

Patch 014
---------
File: ./CMakeLists.txt
Line: 119-120
Reason: Remove compatibility with cuda <5
Comment: lines
    SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-D_CUDA_5PLUS_SDK")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CUDA_5PLUS_SDK")     


Patch 015
---------
File: ./source/tledCUDAHelpers.h
Lines: 68-72
Reason: Remove compatibility with cuda <5
Comment:
    // #ifndef _CUDA_5PLUS_SDK
    // #define tledCheckCUDAErrors(expr) CUDA_SAFE_CALL(expr)
    // #else
    #define tledCheckCUDAErrors(expr) checkCudaErrors(expr)
    // #endif

Patch 016
---------
File: ./source/niftysim.cpp
Lines: 129-138
Reason: Remove compatibility with cuda <5
Comment:
    #ifdef _GPU_
    if (inParams->cudaDeviceId < 0) {
    // #ifndef _CUDA_5PLUS_SDK
    //      inParams->cudaDeviceId = cutGetMaxGflopsDeviceId();
    // #else
    inParams->cudaDeviceId = gpuGetMaxGflopsDeviceId();
    // #endif
    }
    cudaSetDevice(inParams->cudaDeviceId);   
    #endif

Patch 017
---------
File: ./source/tledCUDA_operators.h
Lines: 4-6
Reason: Remove compatibility with cuda <5
Fix: Comment flag _CUDA_5PLUS_SDK  
TODO: Remove

Patch 018
---------
File: ./source/tledCUDA_operators.h
Lines: Multiple lines
Reason: Remove compatibility with cuda <5
Fix: Comment blocks "#ifndef _CUDA_5PLUS_SDK .. #endif" in lines: 20-25, 36-39, 69-74
TODO: Remove

Patch 019
---------
File: ./source/tledCUDA_operators.cu
Lines: 20-22
Reason: Remove compatibility with cuda <5
Fix: Comment flag _CUDA_5PLUS_SDK  
TODO: Remove

Patch 020
---------
File: ./source/tledCUDA_operators.h
Lines: Multiple lines
Reason: Remove compatibility with cuda <5
Fix: Comment blocks "#ifndef _CUDA_5PLUS_SDK .. #endif" in lines: 25-30, 40-50, 62-82, 85-98, 112-127, 144-149, 157-172, 181-194, 197-208, 211-233, 247-252, 267-279, 302-316, 330-334
TODO: Remove

Patch 021
---------
File: ./testing/helpers/tledCUDAUnitTest.h
Lines: 22-26
Reason: Remove compatibility with cuda <5
Fix: Comment flag _CUDA_5PLUS_SDK  
    // #ifndef _CUDA_5PLUS_SDK
    //    cudaSetDevice(cutGetMaxGflopsDeviceId());
    // #else
        cudaSetDevice(gpuGetMaxGflopsDeviceId());
    //#endif
TODO: Remove


Patch 022
---------
File: ./CMakeLists.txt
From Line: 120 to line: 145
Reason: Update CMake file to modern CUDA

Fix summary:
- Keep CUDA support via:
  - `enable_language(CUDA)`
  - `find_package(CUDA REQUIRED)`
  - include `${CUDA_INCLUDE_DIRS}` and `${CUDA_SDK_COMMON_INCLUDE_DIR}`
- Set target architecture explicitly:
  - `set(CMAKE_CUDA_ARCHITECTURES 86)`
- Keep GPU define propagation:
  - `SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-D_GPU_")`
  - `SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GPU_")`


Patch 023
---------
File: ./CMakeLists.txt
Line: 35
Reason: Introduce explicit switch for modern texture API path.

Add:
    option(TLED_USE_TEXTURE_OBJECTS "Use CUDA texture objects" OFF)

Comment:
- Keeps legacy texture references as default for compatibility.
- Enables progressive migration through a build option.


Patch 024
---------
File: ./source/CMakeLists.txt
Lines: CUDA source list and compile definitions
Reason: Wire texture compatibility wrapper into build system.

Changes:
1) Add compatibility translation unit to CUDA source list:
    tledTextureCompatibility

2) Propagate compile definition to target:
    target_compile_definitions(
        tled
        PRIVATE
        $<$<BOOL:${TLED_USE_TEXTURE_OBJECTS}>:
          TLED_USE_TEXTURE_OBJECTS>
    )


Patch 025
---------
File: ./source/tledTextureCompatibility.h
Reason: Add a compatibility layer between legacy texture references and CUDA texture objects.

What was added:
- `tledTexFetch<T>(...)` for:
  - legacy `texture<T,1,cudaReadModeElementType>`
  - modern `cudaTextureObject_t`
- Helpers:
  - `tled::CreateTextureObject<T>(...)`
  - `tled::DestroyTextureObject(...)`
  - `tledBindTexture(...)`
  - `tledBindTextureSymbol(...)`

Impact:
- Kernels can read through one API (`tledTexFetch`) regardless of the selected
  backend.
- Host code binds texture data through one binding API (`tledBindTexture*`).


Patch 026
---------
File: ./source/tledDeviceTextures.cu
Reason: Support both storage models for texture handles.

Changes:
- Under `#ifndef TLED_USE_TEXTURE_OBJECTS`:
  - Keep legacy `texture<...>` declarations.
- Under `#else`:
  - Define device handles as `__device__ cudaTextureObject_t`.

Impact:
- Global texture symbols remain available to kernels in both build modes.


Patch 027
---------
Files:
- ./source/tledSolverGPU_kernels.cu
- ./source/tledSolverGPU_ROM_kernels.cu
- ./source/tledShellSolver_kernels.cu
- ./source/tledNewmarkTimeStepperGPU_kernels.cu
- ./source/tledCentralDifferenceTimeStepperGPU_kernels.cu
Reason: Remove direct dependency on legacy texture fetch syntax in kernels.

Changes:
- Include `tledTextureCompatibility.h`.
- Replace direct `tex1Dfetch(...)` use with `tledTexFetch<T>(...)`.
- Keep symbol declarations dual-path with `#ifndef TLED_USE_TEXTURE_OBJECTS`.

Impact:
- Same kernel code compiles for both legacy texture references and texture objects.


Patch 028
---------
Files:
- ./source/tledSolverGPU.cu
- ./source/tledSolverGPU_ROM.cu
- ./source/tledShellSolverGPU.cu
- ./source/tledTimeStepperGPU.cu
- ./source/tledNewmarkTimeStepperGPU.cu
- ./source/tledCentralDifferenceTimeStepperGPU.cu
Reason: Modernize host-side texture binding calls.

Changes:
- Replace direct `cudaBindTexture(...)` call sites with:
  - `tledBindTextureSymbol(...)` (and compatibility bind helpers)
- Keep `cudaMemcpyToSymbol(...)` for non-texture constant data.

Impact:
- Binding logic becomes backend-agnostic and supports modern CUDA APIs.


Patch 029
---------
File: ./source/tledSurfaceTopology.tpp
Reason: Auxiliary compatibility cleanup required by the CUDA modernization commit.

Comment:
- Small adjustment included in the same migration set; no API or behavior change
  intended at solver level.


Patch 030
---------
Scope: Migration checkpoint closure
Reason: Complete documentation for the commit:
`ad5ddd7 - ubuntu 22.04, vtk 9.1, cuda 12.2`

This patch closes the gap left after Patch 21 and documents the previously
undocumented texture-object compatibility migration centered on:
- `tledTextureCompatibility.h`
- dual texture declarations in device symbols
- wrapper-based binding/fetch in host and kernels
