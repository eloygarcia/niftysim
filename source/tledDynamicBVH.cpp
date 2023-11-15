// =========================================================================
// File:       tledDynamicBVH.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    September 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledDynamicBVH.h"
#include "tledDynamicBVHUpdater.h"
#include "tledRigidMotionBVHUpdaterCPU.h"
#include "tledStaticBVHCreator.h"
#include "tledMovingRigidContactSurfaceCPU.h"
#include "tledBVHXMLImporter.h"
#include "tledOBB.h"
#include "tledAABB.h"

#ifdef GPU_GP_CONTACT
#include "tledDynamicBVHGPU.h"
#endif

#include <iostream>
#include <cstdlib>

void tledDynamicBVH::SetUpdater(tledDynamicBVHUpdater &r_updater) { 
  mp_Updater = &r_updater; 
  mp_Updater->SetBVH(*this);
}

void tledDynamicBVH::SetBVMaxDisplacement(const float margin) {
  if (margin > GetMargin()) {
    std::cerr << "tledDynamicBVH::SetRealMargin: Real margin has to be < BV margin.\n";
    std::abort();
  }

  m_BVMaxDisplacement = margin; 
}

void tledDynamicBVH::SetMargin(const float bvMargin) {
  tledBVH::SetMargin(bvMargin);
  this->SetBVMaxDisplacement(0.95f*GetMargin());
}

tledDynamicBVH::~tledDynamicBVH() {
  if (mp_Updater != NULL) delete mp_Updater;
}

template <class TBV, class TSurface>
static tledDynamicBVH* _CreateBVH(const tledRigidContactSurface &mesh, const float margin) {
  typedef tledDynamicBVHImpl<TSurface, TBV> __BVH;

  TSurface &r_specSurf = const_cast<TSurface&>(static_cast<const TSurface&>(mesh));
  __BVH *p_bvh;
  tledStaticBVHCreator<__BVH> builder;

  p_bvh = new __BVH(r_specSurf);
  p_bvh->SetMargin(margin);
  p_bvh->SetUpdater(*new tledRigidMotionBVHUpdaterCPU<__BVH>(*p_bvh));    
  p_bvh->Init(builder);

  return p_bvh;
}

template <class TBV>
static tledDynamicBVH* _CreateBVHSwitchSurface(const tledRigidContactSurface &mesh, const float margin) {
  tledDynamicBVH *p_bvh = NULL;

  if (mesh.GetNumberOfFacetVertices() == 3) {
    p_bvh = _CreateBVH<TBV, tledMovingRigidContactSurfaceT3CPU>(mesh, margin);
  } else {
    p_bvh = _CreateBVH<TBV, tledMovingRigidContactSurfaceQ4CPU>(mesh, margin);
  }

  return p_bvh;
}

tledDynamicBVH* tledDynamicBVH::CreateBVH(const tledRigidContactSurface &mesh, const std::string &bvType, const float margin, const bool useGPU) {
  tledDynamicBVH *p_bvh = NULL;

  if (useGPU) {
#ifdef GPU_GP_CONTACT
    p_bvh = tledDynamicBVHGPU::CreateBVH(mesh, bvType, margin);
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
static tledDynamicBVH* _ImportBVH(const tledMovingRigidContactSurfaceCPU &mesh, const XMLNode root) {
  typedef tledDynamicBVHImpl<TSurface, TBV> __BVH;

  TSurface &r_specSurf = const_cast<TSurface&>(static_cast<const TSurface&>(mesh));
  __BVH *p_bvh = new __BVH(r_specSurf);
  tledBVHXMLImporter<__BVH> importer;
  
  importer.SetRootNode(root);
  importer.SetOuputObject(*p_bvh);
  importer.Import();
  p_bvh->SetUpdater(*new tledRigidMotionBVHUpdaterCPU<__BVH>(*p_bvh));    
  p_bvh->GetUpdater().Init();

  return p_bvh;
}

template <class TBV>
static tledDynamicBVH* _ImportRigidBVHSwitchMesh(const tledMovingRigidContactSurfaceCPU &mesh, const XMLNode root) {
  tledDynamicBVH *p_bvh = NULL;  

  if (mesh.GetNumberOfFacetVertices() == 3) {
    p_bvh = _ImportBVH<TBV, tledMovingRigidContactSurfaceT3CPU>(mesh, root);
  } else {
    p_bvh = _ImportBVH<TBV, tledMovingRigidContactSurfaceQ4CPU>(mesh, root);
  }

  return p_bvh;
}

static tledDynamicBVH* _ImportRigidBVHSwitchBV(const tledMovingRigidContactSurfaceCPU &mesh, const XMLNode root) {
  tledDynamicBVH *p_bvh = NULL;
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

tledDynamicBVH* tledDynamicBVH::CreateBVH(const tledRigidContactSurface &mesh, const XMLNode root, const bool useGPU) {
  tledDynamicBVH *p_bvh = NULL;

  if (useGPU) {
#ifdef GPU_GP_CONTACT
    p_bvh = tledDynamicBVHGPU::CreateBVH(mesh, root);
#else
    tledFatalFeatureNotEnabledError;
#endif
  } else {
    p_bvh = _ImportRigidBVHSwitchBV(static_cast<const tledMovingRigidContactSurfaceCPU&>(mesh), root);
  }

  return p_bvh;
}
