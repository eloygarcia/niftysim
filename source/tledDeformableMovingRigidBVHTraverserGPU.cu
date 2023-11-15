// =========================================================================
// File:       tledDeformableMovingRigidBVHTraverserGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    September 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableMovingRigidBVHTraverserGPU_CU
#define tledDeformableMovingRigidBVHTraverserGPU_CU

#include "tledDeformableMovingRigidBVHTraverserGPU.h"
#include "tledDeformableMovingRigidBVHTraverserImplGPU.h"
#include "tledMovingRigidContactSurfaceGPU.h"
#include "tledDynamicBVHGPU.h"
#include "tledDeformableContactSurfaceGPU.h"
#include "tledNarrowConeSelfCollisionBVHGPU.h"

template <class TBV>
tledDeformableMovingRigidBVHTraverserGPU* tledDeformableMovingRigidBVHTraverserGPU::_CreateTraverserForBVType(tledSelfCollisionBVH &r_defBVH, const tledDynamicBVH &rigBVH) {
  typedef tledDynamicBVHImplGPU<tledMovingRigidContactSurfaceT3GPU, TBV> __RigBVH;
  typedef tledNarrowConeSelfCollisionBVHImplGPU<tledDeformableContactSurfaceT3GPU, TBV> __BVH;

  tledDeformableMovingRigidBVHTraverserGPU *p_traverser = NULL;

  assert(typeid(__RigBVH) == typeid(rigBVH));
  if (typeid(r_defBVH) == typeid(__BVH)) {
    p_traverser = new tledDeformableMovingRigidBVHTraverserImplGPU<__BVH, __RigBVH>(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
  } else {
    tledFatalError("Unsupported BVH type");
  }

  return p_traverser;
}

tledDeformableMovingRigidBVHTraverserGPU* tledDeformableMovingRigidBVHTraverserGPU::CreateTraverser(tledSelfCollisionBVH &r_defBVH, const tledDynamicBVH &rigBVH) {
  tledDeformableMovingRigidBVHTraverserGPU *p_traverser = NULL;

  if (r_defBVH.GetBVGeometryID() != rigBVH.GetBVGeometryID() || r_defBVH.GetBVHOrder() != rigBVH.GetBVHOrder()) {
    tledFatalError("Rigid/deformable BVH must be of same type.");
  } else if (r_defBVH.GetBVGeometryID() == tledAABB<2>::BVGeometryID && r_defBVH.GetBVHOrder() == 2) {
    typedef tledDynamicBVHImplGPU<tledMovingRigidContactSurfaceT3GPU, tledAABB<2> > __RigBVH;

    assert(typeid(__RigBVH) == typeid(rigBVH));
    p_traverser = _CreateTraverserForBVType<tledAABB<2> >(r_defBVH, rigBVH);
  } else if (r_defBVH.GetBVGeometryID() == tledAABB<4>::BVGeometryID && r_defBVH.GetBVHOrder() == 4) {
    typedef tledDynamicBVHImplGPU<tledMovingRigidContactSurfaceT3GPU, tledAABB<4> > __RigBVH;

    assert(typeid(__RigBVH) == typeid(rigBVH));
    p_traverser = _CreateTraverserForBVType<tledAABB<4> >(r_defBVH, rigBVH);
  } else if (r_defBVH.GetBVGeometryID() == tledOBB<2>::BVGeometryID && r_defBVH.GetBVHOrder() == 2) {
    typedef tledDynamicBVHImplGPU<tledMovingRigidContactSurfaceT3GPU, tledOBB<2> > __RigBVH;

    assert(typeid(__RigBVH) == typeid(rigBVH));
    p_traverser = _CreateTraverserForBVType<tledOBB<2> >(r_defBVH, rigBVH);
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "BVH type not supported/recognised, is " << typeid(r_defBVH).name());
  }

  return p_traverser;
}

#endif
