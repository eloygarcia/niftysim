// =========================================================================
// File:       tledDeformableDeformableBVHTraverserCPU.cpp
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

#include "tledDeformableContactSurfaceCPU.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"
#include "tledDeformableDeformableBVHTraverserCPU.h"
#include "tledNarrowConeSelfCollisionBVH.h"
#include "tledAABB.h"
#include "tledOBB.h"

#ifdef _USE_BOOST_
#include "tledParallelSelfCollisionBVHTraverserCPU.h"
#endif

template <class TBV>
static tledDeformableDeformableBVHTraverserCPU* _CreateTraverserSwitchMesh(tledSelfCollisionBVH &r_bvh) {
  const tledContactSurface &mesh = r_bvh.GetUnspecifiedMesh();

  tledDeformableDeformableBVHTraverserCPU *p_traverser = NULL;

  if (mesh.GetNumberOfFacetVertices() == 3) {
    if (dynamic_cast<const tledDeformableMembraneContactSurfaceCPU*>(&mesh) != NULL) {
      typedef tledSelfCollisionBVHImpl<tledDeformableMembraneContactSurfaceT3CPU, TBV> __BVH;

#ifdef _USE_BOOST_
      p_traverser = new tledParallelSelfCollisionBVHTraverserCPU<tledDeformableDeformableBVHTraverserImplCPU<__BVH> >(static_cast<__BVH&>(r_bvh));
#else
      p_traverser = new tledDeformableDeformableBVHTraverserImplCPU<__BVH>(static_cast<__BVH&>(r_bvh));
#endif
    } else {
      typedef tledSelfCollisionBVHImpl<tledDeformableContactSurfaceT3CPU, TBV> __BVH;

#ifdef _USE_BOOST_
      p_traverser = new tledParallelSelfCollisionBVHTraverserCPU<tledDeformableDeformableBVHTraverserImplCPU<__BVH> >(static_cast<__BVH&>(r_bvh));
#else
      p_traverser = new tledDeformableDeformableBVHTraverserImplCPU<__BVH>(static_cast<__BVH&>(r_bvh));
#endif
    }
  } else {
    if (dynamic_cast<const tledDeformableMembraneContactSurfaceCPU*>(&mesh) != NULL) {
      typedef tledSelfCollisionBVHImpl<tledDeformableMembraneContactSurfaceQ4CPU, TBV> __BVH;

#ifdef _USE_BOOST_
      p_traverser = new tledParallelSelfCollisionBVHTraverserCPU<tledDeformableDeformableBVHTraverserImplCPU<__BVH> >(static_cast<__BVH&>(r_bvh));
#else
      p_traverser = new tledDeformableDeformableBVHTraverserImplCPU<__BVH>(static_cast<__BVH&>(r_bvh));
#endif
    } else {
      typedef tledSelfCollisionBVHImpl<tledDeformableContactSurfaceQ4CPU, TBV> __BVH;

#ifdef _USE_BOOST_
      p_traverser = new tledParallelSelfCollisionBVHTraverserCPU<tledDeformableDeformableBVHTraverserImplCPU<__BVH> >(static_cast<__BVH&>(r_bvh));
#else
      p_traverser = new tledDeformableDeformableBVHTraverserImplCPU<__BVH>(static_cast<__BVH&>(r_bvh));
#endif
    }
  }

  return p_traverser;
}

static tledDeformableDeformableBVHTraverserCPU* _CreateTraverserSwitchBV(tledSelfCollisionBVH &r_bvh) {
  tledDeformableDeformableBVHTraverserCPU *p_traverser = NULL;

  if (r_bvh.GetBVHOrder() == 2) {
    if (r_bvh.GetBVGeometryID() == tledOBB<2>::BVGeometryID) {
      p_traverser = _CreateTraverserSwitchMesh<tledOBB<2> >(r_bvh);
    } else if (r_bvh.GetBVGeometryID() == tledAABB<2>::BVGeometryID) {
      p_traverser = _CreateTraverserSwitchMesh<tledAABB<2> >(r_bvh);
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "BV geometry type " << r_bvh.GetBVGeometryID() << " is not supported.");
    }
  } else if (r_bvh.GetBVHOrder() == 4) {
    if (r_bvh.GetBVGeometryID() == tledOBB<4>::BVGeometryID) {
      p_traverser = _CreateTraverserSwitchMesh<tledOBB<4> >(r_bvh);
    } else if (r_bvh.GetBVGeometryID() == tledAABB<4>::BVGeometryID) {
      p_traverser = _CreateTraverserSwitchMesh<tledAABB<4> >(r_bvh);
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "BV geometry type " << r_bvh.GetBVGeometryID() << " is not supported.");
    }
  } else {
    tledFatalError("BVH orders > 4 are not supported.");
  }

  return p_traverser;
}

tledDeformableDeformableBVHTraverserCPU* tledDeformableDeformableBVHTraverserCPU::CreateTraverser(tledSelfCollisionBVH &r_bvh) {
  return _CreateTraverserSwitchBV(r_bvh);
}
