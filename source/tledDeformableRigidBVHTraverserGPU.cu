// =========================================================================
// File:       tledDeformableRigidBVHTraverserGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableRigidBVHTraverserGPU_CU
#define tledDeformableRigidBVHTraverserGPU_CU

#include "tledDeformableRigidBVHTraverserGPU.h"
#include "tledDeformableRigidBVHTraverserImplGPU.h"
#include "tledDeformableMovingRigidBVHTraverserGPU.h"
#include "tledSelfCollisionBVHGPU.h"
#include "tledDynamicBVH.h"

tledDeformableRigidBVHTraverserGPU* tledDeformableRigidBVHTraverserGPU::CreateTraverser(tledSelfCollisionBVHGPU &r_bvh, const tledBVH &rigBVH) {
  tledDeformableRigidBVHTraverserGPU *p_traverser = NULL;

  if (r_bvh.GetUnspecifiedMesh().GetNumberOfFacetVertices() == rigBVH.GetUnspecifiedMesh().GetNumberOfFacetVertices()) {
    typedef tledBVHImplGPU<tledBVHImpl<tledRigidContactSurfaceT3GPU, tledAABB<2>, tledStaticBVHGPU> > __RigidBVHT3AABB2;
    typedef tledNarrowConeSelfCollisionBVHImplGPU<tledDeformableContactSurfaceT3GPU, tledAABB<2> > __DeformableBVHT3AABB2;

    if (dynamic_cast<__DeformableBVHT3AABB2*>(&r_bvh) != NULL && dynamic_cast<const __RigidBVHT3AABB2*>(&rigBVH) != NULL) {
      p_traverser = new tledDeformableRigidBVHTraverserImplGPU<__DeformableBVHT3AABB2, __RigidBVHT3AABB2>(static_cast<__DeformableBVHT3AABB2&>(r_bvh), static_cast<const __RigidBVHT3AABB2&>(rigBVH));
    } else {
      tledFatalNotYetImplementedError;
    }
  } else {
    tledFatalError("Can only handle rigid & deformable surfaces with same element type.");
  }

  return p_traverser;
}

#endif
