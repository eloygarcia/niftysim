#-----------------------------------------------------------------------------
## LIBRARIES
#-----------------------------------------------------------------------------

SET(NAME tled)

IF(USE_CUDA)
  SET( CPP_SOURCE_FILES_WITH_HEADERS
      tledHelper
      tledElementH8
      tledForceConstraint
      tledMaterialNHV
      tledMesh
      tledConstraintManager
      tledElementT4ANP
      tledMaterial
      tledMaterialTI
      tledModel
      tledParallel
      tledParallelBVHTraverserCPU
      tledDispConstraint
      tledConstraint
      tledElementT4
      tledMaterialLE
      tledMaterialTIV
      tledSolver
      tledElement
      tledFixConstraint
      tledMaterialNH
      tledTimer
      tledSolverCPU
      tledParallelSolverCPU
      tledParallelShellSolverCPU
      tledSolutionWriter
      tledMaterialAB
      tledMaterialPY
      tledSimulator
      tledGravityConstraint
      tledPressureConstraint
      tledSurfaceConstraint
      tledTractionConstraint
      tledMeshTopology
      tledVectorArithmetic
      tledSurface
      tledVTKMeshLoader
      tledMeshLoader
      tledMSHMeshLoader
      tledSubModelManager
      tledContactSurface
      tledDeformableContactSurface
      tledDeformableContactSurfaceCPU
      tledDeformableMembraneContactSurface
      tledDeformableMembraneContactSurfaceCPU
      tledSelfCollisionBVHCreator
      tledBVHTraverserCPU
      tledDeformableDeformableBVHTraverserCPU
      tledDeformableRigidBVHTraverserCPU
      tledSelfCollisionBVHTraverserCPU
      tledDeformableMovingRigidBVHTraverserCPU
      tledMeshSurface
      tledMembraneMaterialLinear
      tledMembraneMaterialNeoHookean
      tledShellMesh
      tledShellSolver
      tledShellSolverCPU
      tledElementMembraneImpl
      tledElementMembraneSimpleLinear
      tledElementMembraneNonLinear
      tledUnstructuredContactManager
      tledDynamicBVH
      tledStaticBVH
      tledSelfCollisionBVH
      tledContactSolver
      tledContactSolverCPU
      tledDeformableDeformableContactSolver
      tledDeformableRigidContactSolver
      tledDeformableMovingRigidContactSolver
      tledDeformableDeformableContactSolverCPU
      tledDeformableRigidContactSolverCPU
      tledDeformableMovingRigidContactSolverCPU
      tledRigidContactSurface
      tledRigidContactSurfaceCPU
      tledMovingRigidContactSurface
      tledMovingRigidContactSurfaceCPU
      tledMatrixFunctions
      tledTimeStepper
      tledTimeStepperCPU
      tledCentralDifferenceTimeStepperCPU
      tledNewmarkTimeStepperCPU
      tledElementShellBSTP1
      tledElementShellBSTP1Edge
      tledNodeRejector
      tledBoxNodeRejector
      tledSphereNodeRejector
  )

  SET( CUDA_SOURCE_FILES_WITH_HEADERS
      tledCUDA_operators
      tledSolverGPU_ROM
      tledContactCylinder
      tledContactManager
      tledContactUSProbe
      tledContactPlate
      tledSolverGPU
      tledCUDAHelpers
      tledCUDAMemoryBlock
      tledUnstructuredContactManager
  )
 
  SET( CUDA_SOURCE_FILES_WITHOUT_HEADERS
      tledSolverGPU_kernels
      tledDeviceConstantMemoryVariables
      tledDeviceTextures

      # BV specialisations
      tledBV
      tledAABB
      tledOBB
  )

  SET( HEADER_FILES_WITHOUT_SOURCE_FILES
      tledDeviceDeclarations
      tledHelper
      tledXMLExporter
      tledXMLImporter
      tledBasicMeshFileReader
      tledBasicMeshCreator
      tledBasicSurfaceCreator
      tledVolumeSurfaceExtractor
      tledSurfaceXMLExporter
      tledSurfaceTopology
      tledBasicMeshFileWriter
      tledGenericVTKMeshExporter
      tledGenericVTKMeshSource
      tledVTKSurfaceExporter
      tledVTKMeshExporter

      # ODE
      tledCentralDifferenceTimeStepper
      tledNewmarkTimeStepper

      # Contact modelling template classes
      tledBV
      tledAABB
      tledOBB
      tledSelfCollisionBV
      tledBVH
      tledNarrowConeSelfCollisionBVH
      tledBVHCreator
      tledMembraneSelfCollisionBVHCreator
      tledDynamicBVHUpdater
      tledBVHTraverser
      tledBVHXMLExporter
      tledBVHXMLImporter
      tledSelfCollisionBVHXMLExporter
      tledSelfCollisionBVHXMLImporter
      tledBVXMLExporter
      tledBVXMLImporter
      tledAABBXMLExporter
      tledAABBXMLImporter
      tledOBBXMLExporter
      tledOBBXMLImporter
      tledSelfCollisionBVXMLExporter
      tledSelfCollisionBVXMLImporter
      tledSurfaceLoader
      tledContactSurfaceCreator
      tledVTKSurfaceLoader
      tledXMLSurfaceCreator
      tledContactVolumeSurfaceExtractor
      tledContactSurfaceXMLExporter
      tledDeformableContactSurfaceXMLExporter
      tledDeformableMembraneContactSurfaceXMLExporter
      tledContactSurfaceCPU
      tledGeometricSelfCollisionBVHUpdater
      tledGeometricSelfCollisionBVHUpdaterCPU
      tledGreedySelfCollisionBVHUpdater

      tledTimeStepperGPU
      tledCentralDifferenceTimeStepperGPU
      tledNewmarkTimeStepperGPU
      )

    IF (USE_GPU_GP_CONTACT)
      LIST(APPEND HEADER_FILES_WITHOUT_SOURCE_FILES
	tledContactSolverGPU
	tledDeformableDeformableContactSolverGPU
	tledDeformableRigidContactSolverGPU
	tledDeformableMovingRigidContactSolverGPU

	tledBVHTraverserGPU
	tledDeformableDeformableBVHTraverserGPU
	tledDeformableMovingRigidBVHTraverserGPU
	tledDeformableRigidBVHTraverserGPU

	tledSelfCollisionBVHGPU
	tledDynamicBVHGPU
	)

      LIST(APPEND CUDA_SOURCE_FILES_WITH_HEADERS 
	tledContactSurfaceGPU
	tledRigidContactSurfaceGPU
	tledMovingRigidContactSurfaceGPU
	tledDeformableContactSurfaceGPU

	tledBVHGPU
	tledStaticBVHGPU
	tledNarrowConeSelfCollisionBVHGPU
	)
    ENDIF (USE_GPU_GP_CONTACT)

ELSE(USE_CUDA)
  SET( CPP_SOURCE_FILES_WITH_HEADERS
      tledElementH8
      tledHelper
      tledForceConstraint
      tledMaterialNHV
      tledMesh
      tledConstraintManager
      tledElementT4ANP
      tledMaterial
      tledMaterialTI
      tledModel
      tledParallel
      tledParallelBVHTraverserCPU
      tledDispConstraint
      tledConstraint
      tledElementT4
      tledMaterialLE
      tledMaterialTIV
      tledSolver
      tledElement
      tledFixConstraint
      tledMaterialNH
      tledTimer
      tledSolverCPU
      tledParallelSolverCPU
      tledParallelShellSolverCPU
      tledSolutionWriter
      tledContactCylinder
      tledContactManager
      tledContactUSProbe
      tledContactPlate
      tledMaterialAB
      tledMaterialPY
      tledSimulator
      tledVectorArithmetic
      tledMeshTopology
      tledGravityConstraint
      tledPressureConstraint
      tledSurfaceConstraint
      tledTractionConstraint
      tledSurface
      tledVTKMeshLoader
      tledMSHMeshLoader
      tledMeshLoader
      tledSubModelManager
      tledMeshSurface
      tledContactSurface
      tledDeformableContactSurface
      tledDeformableContactSurfaceCPU
      tledDeformableMembraneContactSurface
      tledDeformableMembraneContactSurfaceCPU
      tledSelfCollisionBVHCreator
      tledBVHTraverserCPU
      tledDeformableDeformableBVHTraverserCPU
      tledDeformableRigidBVHTraverserCPU
      tledSelfCollisionBVHTraverserCPU
      tledDeformableMovingRigidBVHTraverserCPU
      tledShellMesh
      tledMembraneMaterialLinear
      tledMembraneMaterialNeoHookean
      tledShellSolver
      tledShellSolverCPU
      tledElementMembraneImpl
      tledElementMembraneSimpleLinear
      tledElementMembraneNonLinear
      tledUnstructuredContactManager
      tledDynamicBVH
      tledStaticBVH
      tledSelfCollisionBVH
      tledContactSolver
      tledContactSolverCPU
      tledDeformableDeformableContactSolver
      tledDeformableRigidContactSolver
      tledDeformableMovingRigidContactSolver
      tledDeformableDeformableContactSolverCPU
      tledDeformableRigidContactSolverCPU
      tledDeformableMovingRigidContactSolverCPU
      tledRigidContactSurface
      tledRigidContactSurfaceCPU
      tledMovingRigidContactSurface
      tledMovingRigidContactSurfaceCPU
      tledMatrixFunctions
      tledTimeStepper
      tledTimeStepperCPU
      tledCentralDifferenceTimeStepperCPU
      tledNewmarkTimeStepperCPU
      tledElementShellBSTP1
      tledElementShellBSTP1Edge
      tledNodeRejector
      tledBoxNodeRejector
      tledSphereNodeRejector
      pause
  )

  SET( HEADER_FILES_WITHOUT_SOURCE_FILES
      tledHelper
      tledXMLExporter
      tledXMLImporter
      tledBasicMeshFileReader
      tledBasicMeshCreator
      tledBasicSurfaceCreator
      tledVolumeSurfaceExtractor
      tledSurfaceXMLExporter
      tledSurfaceTopology
      tledBasicMeshFileWriter
      tledGenericVTKMeshExporter
      tledGenericVTKMeshSource
      tledVTKSurfaceExporter
      tledVTKMeshExporter

      # ODE
      tledCentralDifferenceTimeStepper
      tledNewmarkTimeStepper

      # Contact modelling template classes
      tledBV
      tledAABB
      tledOBB
      tledSelfCollisionBV
      tledBVH
      tledNarrowConeSelfCollisionBVH
      tledBVHCreator
      tledMembraneSelfCollisionBVHCreator
      tledDynamicBVHUpdater
      tledBVHTraverser
      tledBVHXMLExporter
      tledBVHXMLImporter
      tledSelfCollisionBVHXMLExporter
      tledSelfCollisionBVHXMLImporter
      tledBVXMLExporter
      tledBVXMLImporter
      tledAABBXMLExporter
      tledAABBXMLImporter
      tledOBBXMLExporter
      tledOBBXMLImporter
      tledSelfCollisionBVXMLExporter
      tledSelfCollisionBVXMLImporter
      tledSurfaceLoader
      tledContactSurfaceCreator
      tledVTKSurfaceLoader
      tledXMLSurfaceCreator
      tledContactVolumeSurfaceExtractor
      tledContactSurfaceXMLExporter
      tledDeformableContactSurfaceXMLExporter
      tledDeformableMembraneContactSurfaceXMLExporter
      tledContactSurfaceCPU
      tledGeometricSelfCollisionBVHUpdater
      tledGeometricSelfCollisionBVHUpdaterCPU
      tledGreedySelfCollisionBVHUpdater
  )
ENDIF(USE_CUDA)

FOREACH( SOURCE_FILE ${CPP_SOURCE_FILES_WITH_HEADERS} )
   SET( TLED_LIST ${TLED_LIST} ${SOURCE_FILE}.cpp )
   INSTALL( FILES ${SOURCE_FILE}.h DESTINATION include )
ENDFOREACH( SOURCE_FILE )

FOREACH( SOURCE_FILE ${CUDA_SOURCE_FILES_WITH_HEADERS} )
   SET( TLED_LIST ${TLED_LIST} ${SOURCE_FILE}.cu )
   INSTALL( FILES ${SOURCE_FILE}.h DESTINATION include )
ENDFOREACH( SOURCE_FILE )

FOREACH( SOURCE_FILE ${CUDA_SOURCE_FILES_WITHOUT_HEADERS} )
   SET( TLED_LIST ${TLED_LIST} ${SOURCE_FILE}.cu )
ENDFOREACH( SOURCE_FILE )

FOREACH( HEADER_FILE ${HEADER_FILES_WITHOUT_SOURCE_FILES} )
   INSTALL( FILES ${HEADER_FILE}.h DESTINATION include )
ENDFOREACH( HEADER_FILE )
#-----------------------------------------------------------------------------

FILE(GLOB NIFTYSIM_TPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")

SET(NON_EXPORT_TPPS
  "${CMAKE_CURRENT_SOURCE_DIR}/tledGreedySelfCollisionBVHUpdaterGPU.tpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/tledBottomUpBVHUpdaterGPU.tpp"
  )

FILE(GLOB KERNEL_TPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*_kernels.tpp")
LIST(APPEND NON_EXPORT_TPPS
  ${KERNEL_TPP_FILES}
  )

FILE(GLOB CONTACT_IMPL_TPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*ImplGPU.tpp")
LIST(APPEND NON_EXPORT_TPPS
  ${CONTACT_IMPL_TPP_FILES}
  )

FOREACH ( TPP ${NON_EXPORT_TPPS} )
  LIST(REMOVE_ITEM NIFTYSIM_TPP_FILES "${TPP}")
ENDFOREACH ( TPP )
#-----------------------------------------------------------------------------

INSTALL(FILES ${NIFTYSIM_TPP_FILES} DESTINATION include)

IF (USE_CUDA)
  CUDA_ADD_LIBRARY( ${NAME} ${TLED_LIST} )
ELSE (USE_CUDA)
  ADD_LIBRARY( ${NAME} ${TLED_LIST} )
ENDIF (USE_CUDA)

IF (USE_BOOST) 
  TARGET_LINK_LIBRARIES(${NAME} xmlParser ${Boost_LIBRARIES})
ELSE (USE_BOOST) 
  TARGET_LINK_LIBRARIES(${NAME} xmlParser)
ENDIF (USE_BOOST) 

INSTALL(TARGETS ${NAME}
   RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
   LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
   ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

#-----------------------------------------------------------------------------

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

#-----------------------------------------------------------------------------
## EXECUTABLES
#-----------------------------------------------------------------------------

IF(USE_CUDA)
   CUDA_ADD_EXECUTABLE(niftysim niftysim.cpp)
   IF(USE_VIZ)
      TARGET_LINK_LIBRARIES(niftysim tled viz)
   ELSE(USE_VIZ)
      TARGET_LINK_LIBRARIES(niftysim tled)
   ENDIF(USE_VIZ)
ELSE(USE_CUDA)
   ADD_EXECUTABLE(niftysim niftysim.cpp)
   IF(USE_VIZ)
      TARGET_LINK_LIBRARIES(niftysim tled viz)
   ELSE(USE_VIZ)
      TARGET_LINK_LIBRARIES(niftysim tled)
   ENDIF(USE_VIZ)
ENDIF(USE_CUDA)
INSTALL_TARGETS(/bin niftysim)

