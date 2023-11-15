// =========================================================================
// File:       tledDeformableMovingRigidBVHTraverserCPU.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledAABB.h"
#include "tledVectorArithmetic.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"
#include "tledNarrowConeSelfCollisionBVH.h"
#include "tledDynamicBVH.h"
#include "tledMovingRigidContactSurfaceCPU.h"
#include "tledDeformableMovingRigidBVHTraverserCPU.h"

#ifdef _USE_BOOST_
#include "tledParallelDeformableMovingRigidBVHTraverserCPU.h"
#endif

#include <typeinfo>

template <class TBV>
tledDeformableMovingRigidBVHTraverserCPU* tledDeformableMovingRigidBVHTraverserCPU::_CreateTraverserForBVType(tledSelfCollisionBVH &r_defBVH, const tledDynamicBVH &rigBVH) {
  typedef tledDynamicBVHImpl<tledMovingRigidContactSurfaceT3CPU, TBV> __RigBVH;

  tledDeformableMovingRigidBVHTraverserCPU *p_traverser = NULL;

  assert(typeid(__RigBVH) == typeid(rigBVH));
  if (typeid(r_defBVH) == typeid(tledNarrowConeSelfCollisionBVH<tledDeformableMembraneContactSurfaceT3CPU, TBV>)) {
    typedef tledNarrowConeSelfCollisionBVH<tledDeformableMembraneContactSurfaceT3CPU, TBV> __BVH;

#ifdef _USE_BOOST_
    p_traverser = new tledParallelDeformableMovingRigidBVHTraverserCPU<tledDeformableMovingRigidBVHTraverserImplCPU<__BVH, __RigBVH> >(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#else
    p_traverser = new tledDeformableMovingRigidBVHTraverserImplCPU<__BVH, __RigBVH>(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#endif
  } else if (typeid(r_defBVH) == typeid(tledNarrowConeSelfCollisionBVH<tledDeformableContactSurfaceT3CPU, TBV>)) {
    typedef tledNarrowConeSelfCollisionBVH<tledDeformableContactSurfaceT3CPU, TBV> __BVH;

#ifdef _USE_BOOST_
    p_traverser = new tledParallelDeformableMovingRigidBVHTraverserCPU<tledDeformableMovingRigidBVHTraverserImplCPU<__BVH, __RigBVH> >(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#else
    p_traverser = new tledDeformableMovingRigidBVHTraverserImplCPU<__BVH, __RigBVH>(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#endif
  } else {
    tledFatalError("Unsupported BVH type");
  }

  return p_traverser;
}

tledDeformableMovingRigidBVHTraverserCPU* tledDeformableMovingRigidBVHTraverserCPU::CreateTraverser(tledSelfCollisionBVH &r_defBVH, const tledDynamicBVH &rigBVH) {
  tledDeformableMovingRigidBVHTraverserCPU *p_traverser = NULL;

  if (r_defBVH.GetBVGeometryID() != rigBVH.GetBVGeometryID() || r_defBVH.GetBVHOrder() != rigBVH.GetBVHOrder()) {
    tledFatalError("Rigid/deformable BVH must be of same type.");
  } else if (r_defBVH.GetBVGeometryID() == tledAABB<2>::BVGeometryID && r_defBVH.GetBVHOrder() == 2) {
    typedef tledDynamicBVHImpl<tledMovingRigidContactSurfaceT3CPU, tledAABB<2> > __RigBVH;

    assert(typeid(__RigBVH) == typeid(rigBVH));
    p_traverser = _CreateTraverserForBVType<tledAABB<2> >(r_defBVH, rigBVH);
  } else if (r_defBVH.GetBVGeometryID() == tledAABB<4>::BVGeometryID && r_defBVH.GetBVHOrder() == 4) {
    typedef tledDynamicBVHImpl<tledMovingRigidContactSurfaceT3CPU, tledAABB<4> > __RigBVH;

    assert(typeid(__RigBVH) == typeid(rigBVH));
    p_traverser = _CreateTraverserForBVType<tledAABB<4> >(r_defBVH, rigBVH);
  } else if (r_defBVH.GetBVGeometryID() == tledOBB<2>::BVGeometryID && r_defBVH.GetBVHOrder() == 2) {
    typedef tledDynamicBVHImpl<tledMovingRigidContactSurfaceT3CPU, tledOBB<2> > __RigBVH;

    assert(typeid(__RigBVH) == typeid(rigBVH));
    p_traverser = _CreateTraverserForBVType<tledOBB<2> >(r_defBVH, rigBVH);
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "BVH type not supported/recognised, is " << typeid(r_defBVH).name());
  }

  return p_traverser;
}
