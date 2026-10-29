## CPP compatibility: From GCC11 to GCC17

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
Issue: Niftysim test are compiled before compiling the main application
Fix: Realocate IF (BUILD_NIFTYSIM_TESTS) .. ENDIF()

Patch 004
---------
File: ./CMakeLists.txt
From Line: 110 to line: 128
Reason: Revome compatibility with CUDA<5.0 in IF(USE_CUDA) .. ENDIF()

Patch 005
---------
File: ./CMakeLists.txt
From Line: 138 to line: 147
Reason: Revome Error Message CUDA<4.2 in IF(USE_GPU_GP_CONTACT) .. ENDIF()

Patch 006
---------
File: ./CMakeLists.txt
From Line:  56 to line: 66
Reason: Updating CMake files to make them compatible with VTK9
TODO: Upadate $VTK_DIR

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
TODO: Upadate $VTK_DIR

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
TODO: Upadate $VTK_DIR

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

Current Code:
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
Fix: Replace GetTupleValue by GetTuple

Old code:
    p_array->GetTupleValue(vInd, tmpTuple);

New code:
    p_array->GetTypedTuple(vInd, tmpTuple);

Patch 010
---------
File: ./CMakeLists.txt
Line: 66 (Within IF(USE_VIZ) .. ENDIF())
Reason: Updating $VTK_DIR in CMake files to make them compatible with VTK9
Fix: Include:

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
Reason: Updating $VTK_DIR in CMake files to make them compatible with VTK9
Fix: Modifying:

    target_include_directories(unitTests
         PUBLIC
         /usr/include/vtk-9.1
    )

by:

    target_include_directories(unitTests
         PUBLIC
         ${NIFTYSIM_VTK_INCLUDE_DIR}
    )

Patch 12
---------
File: ./source/CMakeLists.txt
Line: 454
Reason: Updating $VTK_DIR in CMake files to make them compatible with VTK9
Fix: Modifying:

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

Patch 13
---------
File: ./CMakeLists.txt
Line: 47
Reason: Remove warning incompatiblities between CUDA and VTK
Old code:
    IF (USE_CUDA AND USE_VIZ)
        # The standard FindCUDA.cmake shipping with all versions up to and including 3.1 is broken when used in conjunction with VTK >= 6.0
        SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
    ENDIF (USE_CUDA AND USE_VIZ)