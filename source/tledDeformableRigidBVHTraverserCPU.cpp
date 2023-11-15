// =========================================================================
// File:       tledDeformableRigidBVHTraverserCPU.cpp
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

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include "tledDeformableRigidBVHTraverserCPU.h"
#include "tledVectorArithmetic.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"
#include "tledNarrowConeSelfCollisionBVH.h"
#include "tledMovingRigidContactSurfaceCPU.h"
#include "tledDeformableMovingRigidBVHTraverserCPU.h"
#include "tledAABB.h"
#include "tledOBB.h"
#include "tledHelper.h"
#include "tledStaticBVH.h"

#ifdef _USE_BOOST_
#include "tledParallelDeformableRigidBVHTraverserCPU.h"
#endif

#include <cmath>
#include <cassert>
#include <typeinfo>

template <class TBV>
tledDeformableRigidBVHTraverserCPU* tledDeformableRigidBVHTraverserCPU::_CreateTraverserForBVType(tledSelfCollisionBVH &r_defBVH, const tledBVH &rigBVH) {
  tledDeformableRigidBVHTraverserCPU *p_traverser = NULL;

  if (r_defBVH.GetUnspecifiedMesh().GetNumberOfFacetVertices() == 3) {
    typedef tledBVHImpl<tledRigidContactSurfaceT3CPU, TBV, tledStaticBVH> __RigBVH;

    assert(typeid(__RigBVH) == typeid(rigBVH));
    if (typeid(r_defBVH) == typeid(tledNarrowConeSelfCollisionBVH<tledDeformableMembraneContactSurfaceT3CPU, TBV>)) {
      typedef tledNarrowConeSelfCollisionBVH<tledDeformableMembraneContactSurfaceT3CPU, TBV> __BVH;

#ifdef _USE_BOOST_
      p_traverser = new tledParallelDeformableRigidBVHTraverserCPU<tledDeformableRigidBVHTraverserImplCPU<__BVH, __RigBVH> >(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#else
      p_traverser = new tledDeformableRigidBVHTraverserImplCPU<__BVH, __RigBVH>(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#endif
    } else if (typeid(r_defBVH) == typeid(tledNarrowConeSelfCollisionBVH<tledDeformableContactSurfaceT3CPU, TBV>)) {
      typedef tledNarrowConeSelfCollisionBVH<tledDeformableContactSurfaceT3CPU, TBV> __BVH;

#ifdef _USE_BOOST_
      p_traverser = new tledParallelDeformableRigidBVHTraverserCPU<tledDeformableRigidBVHTraverserImplCPU<__BVH, __RigBVH> >(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#else
      p_traverser = new tledDeformableRigidBVHTraverserImplCPU<__BVH, __RigBVH>(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#endif
    }
  } else {
    typedef tledBVHImpl<tledRigidContactSurfaceQ4CPU, TBV, tledStaticBVH> __RigBVH;

    assert(typeid(__RigBVH) == typeid(rigBVH));
    if (typeid(r_defBVH) == typeid(tledNarrowConeSelfCollisionBVH<tledDeformableMembraneContactSurfaceQ4CPU, TBV>)) {
      typedef tledNarrowConeSelfCollisionBVH<tledDeformableMembraneContactSurfaceQ4CPU, TBV> __BVH;

#ifdef _USE_BOOST_
      p_traverser = new tledParallelDeformableRigidBVHTraverserCPU<tledDeformableRigidBVHTraverserImplCPU<__BVH, __RigBVH> >(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#else
      p_traverser = new tledDeformableRigidBVHTraverserImplCPU<__BVH, __RigBVH>(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#endif
    } else if (typeid(r_defBVH) == typeid(tledNarrowConeSelfCollisionBVH<tledDeformableContactSurfaceQ4CPU, TBV>)) {
      typedef tledNarrowConeSelfCollisionBVH<tledDeformableContactSurfaceQ4CPU, TBV> __BVH;

#ifdef _USE_BOOST_
      p_traverser = new tledParallelDeformableRigidBVHTraverserCPU<tledDeformableRigidBVHTraverserImplCPU<__BVH, __RigBVH> >(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#else
      p_traverser = new tledDeformableRigidBVHTraverserImplCPU<__BVH, __RigBVH>(static_cast<__BVH&>(r_defBVH), static_cast<const __RigBVH&>(rigBVH));
#endif
    }
  }

  return p_traverser;
}

tledDeformableRigidBVHTraverserCPU* tledDeformableRigidBVHTraverserCPU::CreateTraverser(tledSelfCollisionBVH &r_defBVH, const tledBVH &rigBVH) {
  tledDeformableRigidBVHTraverserCPU *p_traverser = NULL;

  if (r_defBVH.GetBVGeometryID() != rigBVH.GetBVGeometryID() || r_defBVH.GetBVHOrder() != rigBVH.GetBVHOrder()) {
    tledFatalError("Rigid/deformable BVH must be of same type.");
  } else if (r_defBVH.GetBVGeometryID() == tledAABB<2>::BVGeometryID && r_defBVH.GetBVHOrder() == 2) {
    p_traverser = _CreateTraverserForBVType<tledAABB<2> >(r_defBVH, rigBVH);
  } else if (r_defBVH.GetBVGeometryID() == tledOBB<2>::BVGeometryID && r_defBVH.GetBVHOrder() == 2) {
    p_traverser = _CreateTraverserForBVType<tledOBB<2> >(r_defBVH, rigBVH);
  } else if (r_defBVH.GetBVGeometryID() == tledAABB<4>::BVGeometryID && r_defBVH.GetBVHOrder() == 4) {
    p_traverser = _CreateTraverserForBVType<tledAABB<4> >(r_defBVH, rigBVH);
  } else if (r_defBVH.GetBVGeometryID() == tledOBB<4>::BVGeometryID && r_defBVH.GetBVHOrder() == 4) {
    p_traverser = _CreateTraverserForBVType<tledOBB<4> >(r_defBVH, rigBVH);
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "BVH type not supported/recognised, is " << typeid(r_defBVH).name());
  }

  return p_traverser;
}
