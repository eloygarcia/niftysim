// =========================================================================
// File:       tledStaticBVH.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledStaticBVH.h"
#include "tledAABB.h"
#include "tledOBB.h"
#include "tledBVHXMLImporter.h"
#include "tledRigidContactSurfaceCPU.h"
#include "tledStaticBVHCreator.h"
#if defined _GPU_ && defined GPU_GP_CONTACT
#include "tledRigidContactSurfaceGPU.h"
#include "tledStaticBVHGPU.h"
#endif

#include <typeinfo>

template <class TBV, class TSurface>
static tledStaticBVH* _CreateBVH(const tledRigidContactSurface &mesh, const float margin) {
  typedef tledBVHImpl<TSurface, TBV, tledStaticBVH> __BVH;

  const TSurface &specSurf = static_cast<const TSurface&>(mesh);

  tledStaticBVH *p_bvh;

  tledStaticBVHCreator<__BVH> builder;

  p_bvh = new __BVH(specSurf);
  p_bvh->SetMargin(margin);
  p_bvh->Init(builder);

  return p_bvh;
}

template <class TBV>
static tledStaticBVH* _CreateBVHSwitchSurface(const tledRigidContactSurface &mesh, const float margin) {
  tledStaticBVH *p_bvh = NULL;

  if (mesh.GetNumberOfFacetVertices() == 3) {
    p_bvh = _CreateBVH<TBV, tledRigidContactSurfaceT3CPU>(mesh, margin);
  } else if (mesh.GetNumberOfFacetVertices() == 4) {
    p_bvh = _CreateBVH<TBV, tledRigidContactSurfaceQ4CPU>(mesh, margin);
  } else {
    tledFatalError("Contact surface not supported.");
  }

  return p_bvh;
}

tledStaticBVH* tledStaticBVH::CreateBVH(const tledRigidContactSurface &mesh, const std::string &bvType, const float margin, const bool useGPU) {
  tledStaticBVH *p_bvh = NULL;

  assert(!mesh.IsMoving());
  if (useGPU) {
#if defined _GPU_ && defined GPU_GP_CONTACT
    if (dynamic_cast<const tledRigidContactSurfaceGPU*>(&mesh) == NULL) {
      tledLogErrorStream(tledHelper::FatalError() << "Use-GPU flag set, but contact surface is not GPU enabled. Type is " << typeid(mesh).name());
    }
    p_bvh = tledStaticBVHGPU::CreateBVH(static_cast<const tledRigidContactSurfaceGPU&>(mesh), margin);
#else
    tledFatalFeatureNotEnabledError;
#endif
  } else {
    if (bvType == "AABB" || bvType == "AABB2") {
      p_bvh = _CreateBVHSwitchSurface<tledAABB<2> >(mesh, margin);
    } else if (bvType == "OBB" || bvType == "OBB2") {
      p_bvh = _CreateBVHSwitchSurface<tledOBB<2> >(mesh, margin);
    } else if (bvType == "AABB4") {
      p_bvh = _CreateBVHSwitchSurface<tledAABB<4> >(mesh, margin);
    } else if (bvType == "OBB4") {
      p_bvh = _CreateBVHSwitchSurface<tledOBB<4> >(mesh, margin);
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "BV type " << bvType << " not recognised.");
    }
  }

  return p_bvh;
}

template <class TBV, class TSurface>
static tledStaticBVH* _ImportBVH(const tledRigidContactSurface &mesh, const XMLNode root) {
  typedef tledBVHImpl<TSurface, TBV, tledStaticBVH> __BVH;

  const TSurface &specSurf = static_cast<const TSurface&>(mesh);

  __BVH *p_bvh = new __BVH(specSurf);
  tledBVHXMLImporter<__BVH> importer;
  
  importer.SetRootNode(root);
  importer.SetOuputObject(*p_bvh);
  importer.Import();

  return p_bvh;
}

template <class TBV>
static tledStaticBVH* _ImportRigidBVHSwitchMesh(const tledRigidContactSurface &mesh, const XMLNode root) {
  tledStaticBVH *p_bvh = NULL;  

  if (mesh.GetNumberOfFacetVertices() == 3) {
    p_bvh = _ImportBVH<TBV, tledRigidContactSurfaceT3CPU>(mesh, root);
  } else {
    p_bvh = _ImportBVH<TBV, tledRigidContactSurfaceQ4CPU>(mesh, root);
  }

  return p_bvh;
}

static tledStaticBVH* _ImportRigidBVHSwitchBV(const tledRigidContactSurface &mesh, const XMLNode root) {
  tledStaticBVH *p_bvh = NULL;
  std::string bvType = "AABB";

  if (root.nChildNode("BVType") > 0) {
    bvType = root.getChildNode("BVType").getText();
  }

  if (bvType == "AABB" || bvType == "AABB2") {
    p_bvh = _ImportRigidBVHSwitchMesh<tledAABB<2> >(mesh, root);
  } else if (bvType == "OBB" || bvType == "OBB2") {
    p_bvh = _ImportRigidBVHSwitchMesh<tledOBB<2> >(mesh, root);
  } else if (bvType == "AABB4") {
    p_bvh = _ImportRigidBVHSwitchMesh<tledAABB<4> >(mesh, root);
  } else if (bvType == "OBB4") {
    p_bvh = _ImportRigidBVHSwitchMesh<tledOBB<4> >(mesh, root);
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "BV type " << bvType << " not supported.");
  }

  return p_bvh;
}

tledStaticBVH* tledStaticBVH::CreateBVH(const tledRigidContactSurface &mesh, const XMLNode root, const bool useGPU) {
  tledStaticBVH *p_bvh = NULL;

  assert(!mesh.IsMoving());
  if (useGPU) {
#if defined _GPU_ && defined GPU_GP_CONTACT
    p_bvh = tledStaticBVHGPU::CreateBVH(static_cast<const tledRigidContactSurfaceGPU&>(mesh), root);
#else
    tledFatalFeatureNotEnabledError;
#endif
  } else {
    p_bvh = _ImportRigidBVHSwitchBV(mesh, root);
  }

  return p_bvh;
}
