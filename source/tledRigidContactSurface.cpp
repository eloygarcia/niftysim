// =========================================================================
// File:       tledRigidContactSurface.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledRigidContactSurface.h"
#include "tledRigidContactSurfaceCPU.h"
#include "tledModel.h"

#ifdef GPU_GP_CONTACT
#include "tledRigidContactSurfaceGPU.h"
#endif

float* tledRigidContactSurface::GetFinalDisplacement(float *p_dst) const {
  std::fill(p_dst, p_dst + 3*this->GetNumberOfNodes(), 0.f);
  
  return p_dst;
}

void tledRigidContactSurface::SetSlaveNodeIndices(const std::vector<int> &slaveNodeIndices, const int numDeformableNodes) {
  m_SlaveMask = std::vector<unsigned char>(numDeformableNodes, static_cast<unsigned char>(false));
  for (std::vector<int>::const_iterator ic_n = slaveNodeIndices.begin(); ic_n < slaveNodeIndices.end(); ic_n++) {
    assert(*ic_n < int(m_SlaveMask.size()) && *ic_n >= 0);
    m_SlaveMask[*ic_n] = static_cast<unsigned char>(true);
  }
}

bool tledRigidContactSurface::IsVTKSurfaceXMLRepresentation(const XMLNode rootNode) {
  return rootNode.nChildNode("VTKSurface") > 0;
}

bool tledRigidContactSurface::IsMovingSurfaceXMLRepresentation(const XMLNode rootNode) {
  return rootNode.nChildNode("Motion") > 0;
}

std::string tledRigidContactSurface::ExtractSurfaceType(const XMLNode meshSpec) {
  XMLNode elTypeSpecNode;
  char const *type = NULL;

  if (meshSpec.nChildNode("Elements") > 0) elTypeSpecNode = meshSpec.getChildNode("Elements");
  else if (meshSpec.nChildNode("VTKSurface") > 0) elTypeSpecNode = meshSpec.getChildNode("VTKSurface");

  if (!elTypeSpecNode.isEmpty()) {
    type = elTypeSpecNode.getAttribute("Type");
  }
   
  if (type == NULL) {
    tledFatalError("Could not find an attribute specifying the element type.");
  }  

  return std::string(type);
}

tledRigidContactSurface* tledRigidContactSurface::CreateSurface(const XMLNode &meshSpec, const bool useGPU) {
  tledRigidContactSurface *p_surf = NULL;

  if (useGPU) {
#ifdef GPU_GP_CONTACT
    p_surf = tledRigidContactSurfaceGPU::CreateSurface(meshSpec);
#else
    tledFatalFeatureNotEnabledError;
#endif
  } else {
    p_surf = tledRigidContactSurfaceCPU::CreateSurface(meshSpec);
  }

  return p_surf;
}

