// =========================================================================
// File:       tledSelfCollisionBVHGPU.cu
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
#ifndef tledSelfCollisionBVHGPU_CU
#define tledSelfCollisionBVHGPU_CU

#include "tledSelfCollisionBVHGPU.h"
#include "tledDeformableMembraneContactSurface.h"
#include "tledDeformableContactSurfaceGPU.h"
#include "tledGeometricSelfCollisionBVHUpdaterDevelGPU.h"
#include "tledNarrowConeSelfCollisionBVHGPU.h"
#include "tledGreedySelfCollisionBVHUpdaterGPU.h"
#include "tledSelfCollisionBVHXMLImporter.h"
#include "tledSelfCollisionBVHCreator.h"
#include "tledAABB.h"
#include "tledOBB.h"

#include <cstdlib>

template <class TBV>
tledSelfCollisionBVHGPU* tledSelfCollisionBVHGPU::_CreateBVHForBVType(tledDeformableContactSurface &r_mesh, const float margin, const float maxDisplacement) {
  tledSelfCollisionBVHGPU *p_bvh = NULL;
  tledDynamicBVHUpdater *p_updater = NULL;
  tledBVHCreator *p_builder = NULL;
  
  if (dynamic_cast<tledDeformableMembraneContactSurface*>(&r_mesh) != NULL) {
    tledFatalNotYetImplementedError;
    // switch (r_mesh.GetNumberOfFacetVertices()) {
    // case 3: {
    //   typedef tledDeformableMembraneContactSurfaceImpl<3> __Mesh;
    //   typedef tledNarrowConeSelfCollisionBVHImplGPU<__Mesh, __BV> __BVH;

    //   p_bvh = new tledNarrowConeSelfCollisionBVHImplGPU<__Mesh, __BV>(static_cast<__Mesh&>(r_mesh));
    //   p_builder = new tledMembraneSelfCollisionBVHCreator<__BVH>();
    //   p_updater = new tledGeometricSelfCollisionBVHUpdaterDevelGPU<__BVH>();

    //   break;
    // }

    // case 4:
    //   std::abort();
    //   break;

    // default:
    //   std::abort();
    // }
  } else {
    if (r_mesh.GetNumberOfFacetVertices() == 3) {
      typedef tledDeformableContactSurfaceT3GPU __Mesh;
      typedef tledNarrowConeSelfCollisionBVHImplGPU<__Mesh, TBV> __BVH;

      p_bvh = new tledNarrowConeSelfCollisionBVHImplGPU<__Mesh, TBV>(static_cast<__Mesh&>(r_mesh));
      p_builder = new tledSelfCollisionBVHCreator<__BVH>();
      //p_updater = new tledGeometricSelfCollisionBVHUpdaterDevelGPU<__BVH>();
      p_updater = new tledGreedySelfCollisionBVHUpdaterGPU<__BVH>();
    } else {
      tledFatalNotYetImplementedError;
    }
  }

  p_bvh->SetMargin(margin);
  p_bvh->SetBVMaxDisplacement(maxDisplacement);
  p_bvh->SetUpdater(*p_updater);
  p_bvh->Init(*p_builder);

  delete p_builder;

  return p_bvh;
}

tledSelfCollisionBVHGPU* tledSelfCollisionBVHGPU::CreateBVH(tledDeformableContactSurface &r_mesh, const std::string &bvType, const float margin, const float maxDisplacement) {
  tledSelfCollisionBVHGPU *p_bvh = NULL;

  if (bvType == "AABB" || bvType == "AABB2") {
    typedef tledAABB<2> __BV;

    p_bvh = _CreateBVHForBVType<__BV>(r_mesh, margin, maxDisplacement);
  } else if (bvType == "AABB4") {
    typedef tledAABB<4> __BV;

    p_bvh = _CreateBVHForBVType<__BV>(r_mesh, margin, maxDisplacement);    
  } else if (bvType == "OBB" || bvType == "OBB2") {
    typedef tledOBB<2> __BV;

    p_bvh = _CreateBVHForBVType<__BV>(r_mesh, margin, maxDisplacement);    
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "BV type " << bvType << " is not yet supported on the GPU");
  }

  return p_bvh;
}

template <class TSurface, class TBV>
static tledSelfCollisionBVHGPU* _ImportBVHWithType(tledDeformableContactSurface &r_mesh, const XMLNode root) {
  typedef tledNarrowConeSelfCollisionBVHImplGPU<TSurface, TBV> __BVH;

  __BVH *p_bvh = new __BVH(static_cast<TSurface&>(r_mesh));
  tledSelfCollisionBVHXMLImporter<__BVH> importer;

  importer.SetRootNode(root);
  importer.SetOuputObject(static_cast<__BVH&>(*p_bvh));
  importer.Import();
  p_bvh->SetUpdater(*new tledGreedySelfCollisionBVHUpdaterGPU<__BVH>());

  return p_bvh;
}

tledSelfCollisionBVHGPU* tledSelfCollisionBVHGPU::CreateBVH(tledDeformableContactSurface &r_mesh, const XMLNode root) {
  tledSelfCollisionBVHGPU *p_bvh;
  std::string bvType = "AABB";

  if (root.nChildNode("BVType") > 0) {
    bvType = root.getChildNode("BVType").getText();
  }

  if (dynamic_cast<tledDeformableContactSurfaceT3GPU*>(&r_mesh) != NULL) {
    typedef tledDeformableContactSurfaceT3GPU __Mesh;

    if (bvType == "AABB" || bvType == "AABB2") {
      typedef tledAABB<2> __BV;

      p_bvh = _ImportBVHWithType<__Mesh, __BV>(r_mesh, root);
    } else if (bvType == "AABB4") {
      typedef tledAABB<4> __BV;

      p_bvh = _ImportBVHWithType<__Mesh, __BV>(r_mesh, root);      
    } else if (bvType == "OBB2" || bvType == "OBB") {
      typedef tledOBB<2> __BV;

      p_bvh = _ImportBVHWithType<__Mesh, __BV>(r_mesh, root);
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "BV type " << bvType << " is not yet supported on the GPU");
    }      
  } else {
    tledFatalNotYetImplementedError;
  }

  return p_bvh;
}

#endif
