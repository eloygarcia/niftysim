// =========================================================================
// File:       tledRigidContactSurfaceCPU.cpp
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
#include "tledRigidContactSurfaceCPU.h"
#include "tledMovingRigidContactSurfaceCPU.h"

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

tledRigidContactSurfaceCPU* tledRigidContactSurfaceCPU::CreateSurface(const XMLNode &meshSpec) {
  tledRigidContactSurfaceCPU *p_surf = NULL;

  if (tledRigidContactSurface::ExtractSurfaceType(meshSpec) == "T3") {
    if (tledRigidContactSurface::IsMovingSurfaceXMLRepresentation(meshSpec)) {
      p_surf = _ReadMovingRigidSurface(new tledMovingRigidContactSurfaceT3CPU(), meshSpec);	
    } else {
      p_surf = _ReadBasicRigidSurface(new tledRigidContactSurfaceT3CPU(), meshSpec);	
    }
  } else if (tledRigidContactSurface::ExtractSurfaceType(meshSpec) == "Q4") {
    if (tledRigidContactSurface::IsMovingSurfaceXMLRepresentation(meshSpec)) {
      p_surf = _ReadMovingRigidSurface(new tledMovingRigidContactSurfaceQ4CPU(), meshSpec);	
    } else {
      p_surf = _ReadBasicRigidSurface(new tledRigidContactSurfaceQ4CPU(), meshSpec);	
    }
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Element type " << tledRigidContactSurface::ExtractSurfaceType(meshSpec) << " is (currently) not supported.");
  }

  return p_surf;
}
