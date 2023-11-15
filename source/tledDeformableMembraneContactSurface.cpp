// =========================================================================
// File:       tledDeformableMembraneContactSurface.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledDeformableMembraneContactSurface.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"
#include "tledHelper.h"

tledDeformableContactSurface* tledDeformableMembraneContactSurface::CreateSurface(const std::string &type, const bool useGPU) {
  tledDeformableContactSurface *p_surf = NULL;

  if (useGPU) {
    tledFatalNotYetImplementedError;
  } else {
    p_surf = tledDeformableMembraneContactSurfaceCPU::CreateSurface(type);
  }

  return p_surf;  
}

tledDeformableContactSurface* tledDeformableMembraneContactSurface::CreateSurface(const tledMesh &volumeMesh, const tledSurface &membraneMesh, const bool useGPU) {
  tledDeformableContactSurface *p_surf = NULL;

  if (useGPU) {
    tledFatalNotYetImplementedError;
  } else {
    p_surf = tledDeformableMembraneContactSurfaceCPU::CreateSurface(volumeMesh, membraneMesh);
  }

  return p_surf;
}

tledDeformableContactSurface* tledDeformableMembraneContactSurface::CreateSurface(const XMLNode xmlRep, const bool useGPU) {
  tledDeformableContactSurface *p_surf = NULL;

  if (useGPU) {
    tledFatalNotYetImplementedError;
  } else {
    p_surf = tledDeformableMembraneContactSurfaceCPU::CreateSurface(xmlRep);
  }

  return p_surf;
}
