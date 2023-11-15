// =========================================================================
// File:       tledRigidContactSurfaceGPU.cu
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
#ifndef tledRigidContactSurfaceGPU_CU
#define tledRigidContactSurfaceGPU_CU
#include "tledRigidContactSurfaceGPU.h"
#include "tledHelper.h"
#include "tledCUDAHelpers.h"
#include "tledMovingRigidContactSurfaceGPU.h"

template <class TSurface>
TSurface* _ReadBasicRigidSurface(TSurface *p_surf, const XMLNode &meshSpec) {
  typename TSurface::BasicSurfaceXMLImporter loader;

  loader.SetOuputObject(*p_surf);
  loader.SetRootNode(meshSpec);
  loader.Import();

  return p_surf;
}

template <class TSurface>
TSurface* _ReadMovingRigidSurface(TSurface *p_surf, const XMLNode &meshSpec) {
  typename TSurface::MovingSurfaceXMLImporter loader;

  loader.SetOuputObject(*p_surf);
  loader.SetRootNode(meshSpec);
  loader.Import();

  return p_surf;
}

tledRigidContactSurfaceGPU* tledRigidContactSurfaceGPU::CreateSurface(const XMLNode &meshSpec) {
  tledRigidContactSurfaceGPU *p_surf = NULL;
    
  if (tledRigidContactSurface::IsMovingSurfaceXMLRepresentation(meshSpec)) {
    if (tledRigidContactSurface::ExtractSurfaceType(meshSpec) == "T3") {
      p_surf = _ReadMovingRigidSurface(new tledMovingRigidContactSurfaceT3GPU(), meshSpec);
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "Element type " << tledRigidContactSurface::ExtractSurfaceType(meshSpec) << " is (currently) not supported.");
    }
  } else {
    if (tledRigidContactSurface::ExtractSurfaceType(meshSpec) == "T3") {
      p_surf = _ReadBasicRigidSurface(new tledRigidContactSurfaceT3GPU(), meshSpec);
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "Element type " << tledRigidContactSurface::ExtractSurfaceType(meshSpec) << " is (currently) not supported.");
    }
  }

  return p_surf;
}

#endif
