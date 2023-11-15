CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/helpers/tledTestResources.h.in ${CMAKE_CURRENT_BINARY_DIR}/helpers/tledTestResources.h @ONLY)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/helpers)

ADD_SUBDIRECTORY(helpers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/helpers)

# List of unit tests that require VTK
SET(VTK_ONLY_TESTS
  testMSHMeshLoader
  testVTKMeshLoader
  testVTKMeshExporter
  testRigidContactSurface
  testMovingRigidContactSurface
)

# List of unit tests for which VTK is optional
SET(NO_VTK_TESTS
  testHelpers
  testModel
  testConstraintManager
  testTractionConstraint
  testPressureConstraint
  testSolver
  testMeshTopology
  testSubModelManager
  testVectorArithmetic
  testDeformableContactSurface
  testMatrixFunctions
  testElementMembraneNonLinear
  testElementMembraneSimpleLinear
  testElementShellBSTP1
  testContactSolverCPU
  testDeformableContactSurfaceXMLExporter
  testDeformableMembraneContactSurfaceXMLExporter
  testSelfCollisionBVHXMLExporter
  testXMLDeformableContactSurfaceCreator
  testXMLDeformableMembraneContactSurfaceCreator
  testSelfCollisionBVHXMLImporter
  testBVHTraverserCPU
  testDeformableContactSurfaceCPU
  testSelfCollisionBVHCreator
  testSelfCollisionBVH
  testOBB
)

# List of pure CUDA tests
SET(CUDA_ONLY_TESTS
  testCUDAHelpers
  testCUDAMemoryBlock
  testComplexCUDAHelpers
  testBVHTraverserGPU
  testDeformableContactSurfaceGPU
)

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

IF (USE_VIZ)
  FOREACH (TEST ${VTK_ONLY_TESTS}) 
    IF (USE_CUDA)
      CUDA_ADD_EXECUTABLE(${TEST} ${TEST}.cpp)
    ELSE (USE_CUDA)
      ADD_EXECUTABLE(${TEST} ${TEST}.cpp)
    ENDIF (USE_CUDA)

    TARGET_LINK_LIBRARIES(${TEST} unitTests tled viz)
    STRING(REGEX REPLACE "^test" "" TESTNAME ${TEST})
    ADD_TEST(NAME ${TESTNAME} COMMAND ${TEST})
  ENDFOREACH (TEST) 
ENDIF (USE_VIZ)

IF (USE_CUDA)
  FOREACH (TEST ${CUDA_ONLY_TESTS})
    CUDA_ADD_EXECUTABLE(${TEST} ${TEST}.cu)

    IF (USE_VIZ)
      TARGET_LINK_LIBRARIES(${TEST} unitTests tled viz)
    ELSE (USE_VIZ)
      TARGET_LINK_LIBRARIES(${TEST} unitTests tled)
    ENDIF (USE_VIZ)  
    STRING(REGEX REPLACE "^test" "" TESTNAME ${TEST})
    ADD_TEST(NAME ${TESTNAME} COMMAND ${TEST})    
  ENDFOREACH (TEST)
ENDIF (USE_CUDA)
