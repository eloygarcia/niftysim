// =========================================================================
// File:       tledDynamicBVHGPU.cu
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
#ifndef tledDynamicBVHGPU_CU
#define tledDynamicBVHGPU_CU

#include "tledDynamicBVHGPU.h"
#include "tledMovingRigidContactSurfaceGPU.h"
#include "tledStaticBVHCreator.h"
#include "tledRigidMotionBVHUpdaterGPU.h"
#include "tledBVHXMLImporter.h"

template <class TBV, class TSurface>
static tledDynamicBVHGPU* _CreateBVH(const TSurface &specSurf, const float margin) {
  typedef tledDynamicBVHImplGPU<TSurface, TBV> __BVH;

  __BVH *p_bvh;
  tledStaticBVHCreator<__BVH> builder;

  p_bvh = new __BVH(specSurf);
  p_bvh->SetMargin(margin);
  //  p_bvh->SetUpdater(*new tledRigidMotionBVHUpdaterGPU<__BVH>(*p_bvh));    
  p_bvh->SetUpdater(*new tledGreedySelfCollisionBVHUpdaterGPU<__BVH>(*p_bvh));    
  p_bvh->Init(builder);

  return p_bvh;
}

template <class TSurface>
static tledDynamicBVHGPU* _CreateBVHSwitchBVType(const TSurface &mesh, const std::string &bvType, const float margin) {
  tledDynamicBVHGPU *p_bvh = NULL;

  if (bvType == "AABB" || bvType == "AABB2") {
    p_bvh = _CreateBVH<tledAABB<2> >(mesh, margin);
  } else if (bvType == "OBB" || bvType == "OBB2") {
    p_bvh = _CreateBVH<tledOBB<2> >(mesh, margin);
  } else if (bvType == "AABB4") {
    p_bvh = _CreateBVH<tledAABB<4> >(mesh, margin);
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "BV type " << bvType << " not supported.");
  }

  return p_bvh;
}

template <class TBV, class TSurface>
static tledDynamicBVHGPU* _ImportBVH(const tledMovingRigidContactSurfaceGPU &mesh, const XMLNode root) {
  typedef tledDynamicBVHImplGPU<TSurface, TBV> __BVH;

  const TSurface &specSurf = static_cast<const TSurface&>(mesh);

  __BVH *p_bvh = new __BVH(specSurf);
  tledBVHXMLImporter<__BVH> importer;
  
  importer.SetRootNode(root);
  importer.SetOuputObject(*p_bvh);
  importer.Import();
  p_bvh->SetUpdater(*new tledRigidMotionBVHUpdaterGPU<__BVH>(*p_bvh));    
  p_bvh->GetUpdater().Init();

  return p_bvh;
}

template <class TBV>
static tledDynamicBVHGPU* _ImportRigidBVHSwitchMesh(const tledMovingRigidContactSurfaceGPU &mesh, const XMLNode root) {
  tledDynamicBVHGPU *p_bvh = NULL;  

  if (mesh.GetNumberOfFacetVertices() == 3) {
    p_bvh = _ImportBVH<TBV, tledMovingRigidContactSurfaceT3GPU>(mesh, root);
  } else {
    p_bvh = _ImportBVH<TBV, tledMovingRigidContactSurfaceQ4GPU>(mesh, root);
  }

  return p_bvh;
}

static tledDynamicBVHGPU* _ImportRigidBVHSwitchBV(const tledMovingRigidContactSurfaceGPU &mesh, const XMLNode root) {
  tledDynamicBVHGPU *p_bvh = NULL;
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

tledDynamicBVHGPU* tledDynamicBVHGPU::CreateBVH(const tledRigidContactSurface &mesh, const std::string &bvType, const float margin) {
  tledDynamicBVHGPU *p_bvh = NULL;

  assert(dynamic_cast<const tledMovingRigidContactSurfaceGPU*>(&mesh) != NULL);
  if (mesh.GetNumberOfFacetVertices() == 3) {
    p_bvh = _CreateBVHSwitchBVType(static_cast<const tledMovingRigidContactSurfaceT3GPU&>(mesh), bvType, margin);
  } else {
    tledFatalNotYetImplementedError;
  }

  return p_bvh;
}

tledDynamicBVHGPU* tledDynamicBVHGPU::CreateBVH(const tledRigidContactSurface &mesh, const XMLNode root) {
  assert(dynamic_cast<const tledMovingRigidContactSurfaceGPU*>(&mesh) != NULL);
  return _ImportRigidBVHSwitchBV(static_cast<const tledMovingRigidContactSurfaceGPU&>(mesh), root);
}

#endif



