// =========================================================================
// File:       tledSelfCollisionBVH.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2014
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledSelfCollisionBVH.h"
#include "tledNarrowConeSelfCollisionBVH.h"
#include "tledAABB.h"
#include "tledOBB.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"
#include "tledSelfCollisionBVHXMLImporter.h"
#include "tledGreedySelfCollisionBVHUpdater.h"
#include "tledGeometricSelfCollisionBVHUpdaterCPU.h"
#include "tledMembraneSelfCollisionBVHCreator.h"
#ifdef _USE_BOOST_
#include "tledParallelGeometricSelfCollisionBVHUpdaterCPU.h"
#include "tledParallelSelfCollisionBVHCreator.h"
#endif
#ifdef GPU_GP_CONTACT
#include "tledSelfCollisionBVHGPU.h"
#include "tledDeformableContactSurfaceGPU.h"
#endif

template <class TBVH>
static void _SetSelfCollisionBVHUpdater(TBVH &r_bvh) {
  if (r_bvh.GetMesh().GetNumberOfFacets() > 100) {
#ifdef _USE_BOOST_
    r_bvh.SetUpdater(*new tledParallelGeometricSelfCollisionBVHUpdaterCPU<TBVH>());
#else
    r_bvh.SetUpdater(*new tledGeometricSelfCollisionBVHUpdaterCPU<TBVH>());
#endif
  } else r_bvh.SetUpdater(*new tledGreedySelfCollisionBVHUpdater<TBVH>());
}

template <class TBVH, class TBuilder>
static tledSelfCollisionBVH* _CreateNarrowConeDeformableBVH(typename TBVH::ContactMesh &r_surface, const float margin, const float rateDist) {
  TBVH *p_bvh = new TBVH(r_surface);
  TBuilder builder;

  _SetSelfCollisionBVHUpdater(*p_bvh);

  p_bvh->SetMargin(margin);
  p_bvh->SetBVMaxDisplacement(p_bvh->GetMargin() - 1.5f*rateDist);
  p_bvh->Init(builder);

  return p_bvh;
}

template <class TBV>
static tledSelfCollisionBVH* _CreateNarrowConeDeformableBVHSwitchSurface(tledDeformableContactSurface &r_mesh, const float margin, const float rateDist) {
  tledSelfCollisionBVH *p_bvh = NULL;

  if (r_mesh.GetNumberOfFacetVertices() == 3) {
    if (dynamic_cast<tledDeformableMembraneContactSurfaceCPU*>(&r_mesh) != NULL) {
      typedef tledDeformableMembraneContactSurfaceT3CPU __Mesh;
      typedef tledNarrowConeSelfCollisionBVH<__Mesh, TBV> __BVH;

      p_bvh = _CreateNarrowConeDeformableBVH<__BVH, tledMembraneSelfCollisionBVHCreator<__BVH> >(static_cast<__Mesh&>(r_mesh), margin, rateDist);
    } else {
      typedef tledDeformableContactSurfaceT3CPU __Mesh;
      typedef tledNarrowConeSelfCollisionBVH<__Mesh, TBV> __BVH;
#ifdef _USE_BOOST_
      typedef tledParallelSelfCollisionBVHCreator<__BVH> __Builder;
#else 
      typedef tledSelfCollisionBVHCreator<__BVH> __Builder;
#endif

      p_bvh = _CreateNarrowConeDeformableBVH<__BVH,  __Builder>(static_cast<__Mesh&>(r_mesh), margin, rateDist);
    }
  } else {
    if (dynamic_cast<tledDeformableMembraneContactSurfaceCPU*>(&r_mesh) != NULL) {
      typedef tledDeformableMembraneContactSurfaceQ4CPU __Mesh;
      typedef tledNarrowConeSelfCollisionBVH<__Mesh, TBV> __BVH;

      p_bvh = _CreateNarrowConeDeformableBVH<__BVH, tledMembraneSelfCollisionBVHCreator<__BVH> >(static_cast<__Mesh&>(r_mesh), margin, rateDist);
    } else {
      typedef tledDeformableContactSurfaceQ4CPU __Mesh;
      typedef tledNarrowConeSelfCollisionBVH<__Mesh, TBV> __BVH;
#ifdef _USE_BOOST_
typedef tledParallelSelfCollisionBVHCreator<__BVH> __Builder;
#else 
typedef tledSelfCollisionBVHCreator<__BVH> __Builder;
#endif

      p_bvh = _CreateNarrowConeDeformableBVH<__BVH, __Builder>(static_cast<__Mesh&>(r_mesh), margin, rateDist);
    }
  }

  return p_bvh;
}

tledSelfCollisionBVH* tledSelfCollisionBVH::CreateBVH(tledDeformableContactSurface &r_mesh, const std::string &bvType, const float margin, const float maxDisplacement, const bool useGPU) {
  tledSelfCollisionBVH *p_bvh = NULL;

  if (useGPU) {
#if defined _GPU_ && defined GPU_GP_CONTACT
    p_bvh = tledSelfCollisionBVHGPU::CreateBVH(r_mesh, bvType, margin, maxDisplacement);
#else
    tledFatalFeatureNotEnabledError;
#endif
  } else {
    if (bvType == "AABB" || bvType == "AABB2") {
      p_bvh = _CreateNarrowConeDeformableBVHSwitchSurface<tledAABB<2> >(r_mesh, margin, maxDisplacement);
    } else if (bvType == "OBB" || bvType == "OBB2") {
      p_bvh = _CreateNarrowConeDeformableBVHSwitchSurface<tledOBB<2> >(r_mesh, margin, maxDisplacement);
    } else if (bvType == "AABB4") {
      p_bvh = _CreateNarrowConeDeformableBVHSwitchSurface<tledAABB<4> >(r_mesh, margin, maxDisplacement);
    } else if (bvType == "OBB4") {
      p_bvh = _CreateNarrowConeDeformableBVHSwitchSurface<tledOBB<4> >(r_mesh, margin, maxDisplacement);
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "BV type " << bvType << " not supported.");
    }
  }

  return p_bvh;
}

template <class TBVH>
static tledSelfCollisionBVH* _ImportDeformableBVH(typename TBVH::ContactMesh &r_surface, const XMLNode root) {
  TBVH *p_bvh = new TBVH(r_surface);
  tledSelfCollisionBVHXMLImporter<TBVH> importer;

  importer.SetRootNode(root);
  importer.SetOuputObject(*p_bvh);
  importer.Import();
  _SetSelfCollisionBVHUpdater(*p_bvh);

  return p_bvh;
}

template <class TBV>
static tledSelfCollisionBVH* _ImportDeformableBVHSwitchSurface(tledDeformableContactSurface &r_surface, const XMLNode root) {
  tledSelfCollisionBVH *p_bvh = NULL;

  if (r_surface.GetNumberOfFacetVertices() == 3) {    
    if (dynamic_cast<tledDeformableMembraneContactSurface*>(&r_surface) != NULL) {
      typedef tledNarrowConeSelfCollisionBVH<tledDeformableMembraneContactSurfaceT3CPU, TBV> __BVH;

      p_bvh = _ImportDeformableBVH<__BVH>(static_cast<tledDeformableMembraneContactSurfaceT3CPU&>(r_surface), root);
    } else {
      typedef tledNarrowConeSelfCollisionBVH<tledDeformableContactSurfaceT3CPU, TBV> __BVH;

      p_bvh = _ImportDeformableBVH<__BVH>(static_cast<tledDeformableContactSurfaceT3CPU&>(r_surface), root);	
    }
  } else {
    if (dynamic_cast<tledDeformableMembraneContactSurface*>(&r_surface) != NULL) {
      typedef tledNarrowConeSelfCollisionBVH<tledDeformableMembraneContactSurfaceQ4CPU, TBV> __BVH;

      p_bvh = _ImportDeformableBVH<__BVH>(static_cast<tledDeformableMembraneContactSurfaceQ4CPU&>(r_surface), root);
    } else {
      typedef tledNarrowConeSelfCollisionBVH<tledDeformableContactSurfaceQ4CPU, TBV> __BVH;

      p_bvh = _ImportDeformableBVH<__BVH>(static_cast<tledDeformableContactSurfaceQ4CPU&>(r_surface), root);
    }
  }

  return p_bvh;
}

tledSelfCollisionBVH* tledSelfCollisionBVH::CreateBVH(tledDeformableContactSurface &r_mesh, const XMLNode rootNode, const bool useGPU) {
  tledSelfCollisionBVH *p_bvh = NULL;

  if (useGPU) {
#if defined _GPU_ && defined GPU_GP_CONTACT
    p_bvh = tledSelfCollisionBVHGPU::CreateBVH(static_cast<tledDeformableContactSurfaceGPU&>(r_mesh), rootNode);
#else
    tledFatalFeatureNotEnabledError;
#endif
  } else {
    std::string bvType = "AABB";

    if (rootNode.nChildNode("BVType") > 0) {
      bvType = rootNode.getChildNode("BVType").getText();
    }

    if (bvType == "AABB" || bvType == "AABB2") {
      p_bvh = _ImportDeformableBVHSwitchSurface<tledAABB<2> >(r_mesh, rootNode);
    } else if (bvType == "OBB" || bvType == "OBB2") {
      p_bvh = _ImportDeformableBVHSwitchSurface<tledOBB<2> >(r_mesh, rootNode);
    } else if (bvType == "AABB4") {
      p_bvh = _ImportDeformableBVHSwitchSurface<tledAABB<4> >(r_mesh, rootNode);
    } else if (bvType == "OBB4") {
      p_bvh = _ImportDeformableBVHSwitchSurface<tledOBB<4> >(r_mesh, rootNode);
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "BV type " << bvType << " not supported.");
    }
  }

  return p_bvh;
}
