// =========================================================================
// File:       tledSelfCollisionBVHTraverserCPU.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledSelfCollisionBVHTraverserCPU.h"
#include "tledAABB.h"
#include "tledOBB.h"
#include "tledDeformableContactSurfaceCPU.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"
#include "tledNarrowConeSelfCollisionBVH.h"

#ifdef _USE_BOOST_
#include "tledParallelSelfCollisionBVHTraverserCPU.h"
#endif

#include <typeinfo>

template <class TBV>
static tledSelfCollisionBVHTraverserCPU* _CreateTraverserSwitchMesh(tledSelfCollisionBVH &r_bvh) {
  const tledContactSurface &mesh = r_bvh.GetUnspecifiedMesh();

  tledSelfCollisionBVHTraverserCPU *p_traverser = NULL;
  
  if (mesh.GetNumberOfFacetVertices() == 3) {
    if (dynamic_cast<const tledDeformableMembraneContactSurfaceCPU*>(&mesh) == NULL) {
      typedef tledSelfCollisionBVHImpl<tledDeformableContactSurfaceT3CPU, TBV> __BVH;

#ifdef _USE_BOOST_
      return new tledParallelSelfCollisionBVHTraverserCPU<tledSelfCollisionBVHTraverserImplCPU<__BVH> >(static_cast<__BVH&>(r_bvh));
#else
      return new tledSelfCollisionBVHTraverserImplCPU<__BVH>(static_cast<__BVH&>(r_bvh));
#endif
    } else {
    typedef tledSelfCollisionBVHImpl<tledDeformableMembraneContactSurfaceT3CPU, TBV> __BVH;

#ifdef _USE_BOOST_
    return new tledParallelSelfCollisionBVHTraverserCPU<tledSelfCollisionBVHTraverserImplCPU<__BVH> >(static_cast<__BVH&>(r_bvh));
#else
    return new tledSelfCollisionBVHTraverserImplCPU<__BVH>(static_cast<__BVH&>(r_bvh));
#endif
    }
  } else {
    if (dynamic_cast<const tledDeformableMembraneContactSurfaceCPU*>(&mesh) == NULL) {
      typedef tledSelfCollisionBVHImpl<tledDeformableContactSurfaceQ4CPU, TBV> __BVH;

#ifdef _USE_BOOST_
      return new tledParallelSelfCollisionBVHTraverserCPU<tledSelfCollisionBVHTraverserImplCPU<__BVH> >(static_cast<__BVH&>(r_bvh));
#else
      return new tledSelfCollisionBVHTraverserImplCPU<__BVH>(static_cast<__BVH&>(r_bvh));
#endif
    } else {
      typedef tledSelfCollisionBVHImpl<tledDeformableMembraneContactSurfaceQ4CPU, TBV> __BVH;

#ifdef _USE_BOOST_
      return new tledParallelSelfCollisionBVHTraverserCPU<tledSelfCollisionBVHTraverserImplCPU<__BVH> >(static_cast<__BVH&>(r_bvh));
#else
      return new tledSelfCollisionBVHTraverserImplCPU<__BVH>(static_cast<__BVH&>(r_bvh));
#endif
    }
  }

  return p_traverser;
}

static tledSelfCollisionBVHTraverserCPU* _CreateTraverserSwitchBV(tledSelfCollisionBVH &r_bvh) {
  tledSelfCollisionBVHTraverserCPU *p_traverser = NULL;

  if (r_bvh.GetBVHOrder() == 2) {
    if (r_bvh.GetBVGeometryID() == tledAABB<2>::BVGeometryID) {
      return _CreateTraverserSwitchMesh<tledAABB<2> >(r_bvh);
    } else if (r_bvh.GetBVGeometryID() == tledOBB<2>::BVGeometryID) {
      return _CreateTraverserSwitchMesh<tledOBB<2> >(r_bvh);
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "BV geometry type " << r_bvh.GetBVGeometryID() << " not supported.");
    }
  } else if (r_bvh.GetBVHOrder() == 4) {
    if (r_bvh.GetBVGeometryID() == tledAABB<4>::BVGeometryID) {
      return _CreateTraverserSwitchMesh<tledAABB<4> >(r_bvh);
    } else if (r_bvh.GetBVGeometryID() == tledOBB<4>::BVGeometryID) {
      return _CreateTraverserSwitchMesh<tledOBB<4> >(r_bvh);
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "BV geometry type " << r_bvh.GetBVGeometryID() << " not supported.");
    }
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "BVH order " << r_bvh.GetBVHOrder() << " not supported.");
  }

  return p_traverser;
}

tledSelfCollisionBVHTraverserCPU* tledSelfCollisionBVHTraverserCPU::CreateTraverser(tledSelfCollisionBVH &r_bvh) {
  return _CreateTraverserSwitchBV(r_bvh);
}
