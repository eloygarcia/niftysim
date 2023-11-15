// =========================================================================
// File:       tledUnstructuredContactManager.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    September 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledUnstructuredContactManager_CU
#define tledUnstructuredContactManager_CU

#include "tledUnstructuredContactManager.h"
#include "tledCUDAHelpers.h"

#include "tledCUDA_operators.cu"
#include "tledUnstructuredContactManager_kernels.cu"

#ifdef GPU_GP_CONTACT
#include "tledBottomUpBVHUpdaterGPU.h"
#include "tledGreedySelfCollisionBVHUpdaterGPU.h"
#include "tledDeformableDeformableContactSolverGPU.h"
#include "tledDeformableMovingRigidContactSolverGPU.h"
#include "tledDeformableRigidContactSolverGPU.h"

/* Too many inter-dependencies (inlined) -> must be included, not suitable for stand-alone compilation. */
#include "tledOBB_kernels.cu"
#include "tledAABB_kernels.cu"

#include "tledRigidContactSurfaceGPU_kernels.cu"
#include "tledContactSurfaceGPU_kernels.cu"

#include "tledBVHTraverserGPU.cu"
#include "tledDeformableDeformableBVHTraverserGPU.cu"
#include "tledDeformableRigidBVHTraverserGPU.cu"
#include "tledDeformableMovingRigidBVHTraverserGPU.cu"
#include "tledDeformableDeformableBVHTraverserGPU_kernels.cu"

#include "tledSelfCollisionBVHGPU.cu"
#include "tledDynamicBVHGPU.cu"

#include "tledDeformableDeformableContactSolverGPU.cu"
#include "tledDeformableMovingRigidContactSolverGPU.cu"
#include "tledDeformableRigidContactSolverGPU.cu"
#endif

void tledUnstructuredContactManager::_InitGPUSetup() {
  tledCheckCUDAErrors(cudaMemcpyToSymbol(tledUnstructuredContactManager_kernels::c_SurfaceCloseDistance, &m_CloseDistance, sizeof(float)));
}

#ifdef GPU_GP_CONTACT
bool tledUnstructuredContactManager::ComputeDeformableDeformableContactResponses(float4 *dp_R, const float4 *dpc_UNext, const float4 *dpc_UCurr) { 
  return static_cast<tledDeformableDeformableContactSolverGPU*>(mp_DeformableDeformableSolver)->ComputeContactResponses(dp_R, dpc_UNext, dpc_UCurr); 
}

bool tledUnstructuredContactManager::ComputeDeformableRigidContactResponses(float4 *dp_R, const float4 *dpc_UNext, const float4 *dpc_UCurr) {
  bool hadContacts = false;

  assert(this->UseGPU());
  for (std::vector<tledDeformableRigidContactSolver*>::iterator i_solver = mvp_DeformableRigidContactSolvers.begin(); i_solver < mvp_DeformableRigidContactSolvers.end(); i_solver++) {
    hadContacts |= static_cast<tledDeformableRigidContactSolverGPU*>(*i_solver)->ComputeContactResponses(dp_R, dpc_UNext, dpc_UCurr);
  }

  return hadContacts;
}

#endif

#endif
