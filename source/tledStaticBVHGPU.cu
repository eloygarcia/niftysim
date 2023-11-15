// =========================================================================
// File:       tledStaticBVHGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledStaticBVHGPU_CU
#define tledStaticBVHGPU_CU

#include "tledStaticBVHGPU.h"
#include "tledStaticBVHCreator.h"
#include "tledBVHXMLImporter.h"

tledStaticBVHGPU* tledStaticBVHGPU::CreateBVH(const tledRigidContactSurfaceGPU &mesh, const float margin) {
  tledStaticBVHGPU *p_bvh = NULL;  

  if (dynamic_cast<const tledRigidContactSurfaceT3GPU*>(&mesh) != NULL) {
    typedef tledBVHImplGPU<tledBVHImpl<tledRigidContactSurfaceT3GPU, tledAABB<2>, tledStaticBVHGPU> > __BVH;
    tledStaticBVHCreator<__BVH> builder;

    /* Const cast only performed for compat reasons, no write access performed */
    p_bvh = new __BVH(static_cast<tledRigidContactSurfaceT3GPU&>(const_cast<tledRigidContactSurfaceGPU&>(mesh)));
    p_bvh->SetMargin(margin);
    p_bvh->Init(builder);
  } else {
    tledFatalNotYetImplementedError;
  }

  return p_bvh;
}

tledStaticBVHGPU* tledStaticBVHGPU::CreateBVH(const tledRigidContactSurfaceGPU &mesh, const XMLNode root) {
  tledStaticBVHGPU *p_bvh = NULL;

  if (dynamic_cast<const tledRigidContactSurfaceT3GPU*>(&mesh) != NULL) {
    typedef tledBVHImplGPU<tledBVHImpl<tledRigidContactSurfaceT3GPU, tledAABB<2>, tledStaticBVHGPU> > __BVH;

    tledBVHXMLImporter<__BVH> importer;
    
    p_bvh = new __BVH(static_cast<tledRigidContactSurfaceT3GPU&>(const_cast<tledRigidContactSurfaceGPU&>(mesh)));
    importer.SetRootNode(root);
    importer.SetOuputObject(static_cast<__BVH&>(*p_bvh));
    importer.Import();  
  }  
   
  return p_bvh;
}
#endif
