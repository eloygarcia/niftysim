// =========================================================================
// File:       tledDeformableDeformableBVHTraverserGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableDeformableBVHTraverserGPU_CU
#define tledDeformableDeformableBVHTraverserGPU_CU

#include "tledDeformableDeformableBVHTraverserGPU.h"
#include "tledDeformableDeformableBVHTraverserImplGPU.h"
#include "tledDeformableContactSurfaceGPU.h"
#include "tledNarrowConeSelfCollisionBVHGPU.h"
#include "tledDeformableMembraneContactSurface.h"

tledDeformableDeformableBVHTraverserGPU* tledDeformableDeformableBVHTraverserGPU::CreateTraverser(tledSelfCollisionBVH &r_bvh) {
  tledDeformableDeformableBVHTraverserGPU *p_traverser = NULL;

  if (dynamic_cast<const tledDeformableMembraneContactSurface*>(&r_bvh.GetUnspecifiedMesh()) != NULL) {
    /*    typedef tledDeformableMembraneContactSurfaceImpl<3> __Mesh;
    typedef tledNarrowConeSelfCollisionBVHImplGPU<__Mesh, tledAABB<2> > __BVH;

    p_traverser = new tledDeformableDeformableBVHTraverserImplGPU<__BVH>(static_cast<__BVH&>(r_bvh)); */
    tledFatalNotYetImplementedError;
  } else {
    typedef tledDeformableContactSurfaceT3GPU __Mesh;

    if (r_bvh.GetBVGeometryID() == tledAABB<2>::BVGeometryID) {
      if (r_bvh.GetBVHOrder() == 2) {
	typedef tledNarrowConeSelfCollisionBVHImplGPU<__Mesh, tledAABB<2> > __BVH;

	p_traverser = new tledDeformableDeformableBVHTraverserImplGPU<__BVH>(static_cast<__BVH&>(r_bvh));
      } else if (r_bvh.GetBVHOrder() == 4) {
	typedef tledNarrowConeSelfCollisionBVHImplGPU<__Mesh, tledAABB<4> > __BVH;

	p_traverser = new tledDeformableDeformableBVHTraverserImplGPU<__BVH>(static_cast<__BVH&>(r_bvh));
      } else {
	tledFatalNotYetImplementedError;
      }
    } else {
      tledFatalNotYetImplementedError;
    }
  }

  return p_traverser;
}

#endif
