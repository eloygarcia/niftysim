// =========================================================================
// File:       tledDeformableContactSurface.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledDeformableContactSurface.h"
#include "tledContactVolumeSurfaceExtractor.h"
#include "tledDeformableContactSurfaceCPU.h"
#include "tledVectorArithmetic.h"

#ifdef GPU_GP_CONTACT
#include "tledDeformableContactSurfaceGPU.h"
#endif

tledDeformableContactSurface* tledDeformableContactSurface::CreateSurface(const tledMesh &mesh, const bool useGPU) {
#ifdef GPU_GP_CONTACT
  if (useGPU) {
    return tledDeformableContactSurfaceGPU::CreateSurface(mesh);
  } else {    
#endif
    return tledDeformableContactSurfaceCPU::CreateSurface(mesh);
#ifdef GPU_GP_CONTACT
  }
#endif
}

tledDeformableContactSurface* tledDeformableContactSurface::CreateSurface(const std::string &type, const bool useGPU) {
  tledDeformableContactSurface *p_surf = NULL;

#ifdef GPU_GP_CONTACT
  if (useGPU) {
    p_surf = tledDeformableContactSurfaceGPU::CreateSurface(type);
  } else {
#endif
    p_surf = tledDeformableContactSurfaceCPU::CreateSurface(type);
#ifdef GPU_GP_CONTACT
  }
#endif

  return p_surf;
}

tledDeformableContactSurface* tledDeformableContactSurface::CreateSurface(const XMLNode xmlRep, const bool useGPU) {
  tledDeformableContactSurface *p_surf = NULL;

#ifdef GPU_GP_CONTACT
  if (useGPU) {
    p_surf = tledDeformableContactSurfaceGPU::CreateSurface(xmlRep);
  } else {
#endif
    p_surf = tledDeformableContactSurfaceCPU::CreateSurface(xmlRep);
#ifdef GPU_GP_CONTACT
  }
#endif

  return p_surf;
}

void tledDeformableContactSurface::InitNodeMasses(const float globalMasses[]) {
  m_SurfaceNodeMasses.clear();
  m_SurfaceNodeMasses.reserve(GetSurface2VolumeNodeMap().size());
  for (std::vector<int>::const_iterator ic_ind = GetSurface2VolumeNodeMap().begin(); ic_ind < GetSurface2VolumeNodeMap().end(); ic_ind++) m_SurfaceNodeMasses.push_back(globalMasses[*ic_ind]);
}
