// =========================================================================
// File:       tledMovingRigidContactSurfaceCPU.cpp
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
#include "tledMovingRigidContactSurfaceCPU.h"
#include "tledHelper.h"

#include <algorithm>
#include <cmath>

void tledMovingRigidContactSurfaceCPU::RotateVector(float *p_v, const float rc[], const float rs[]) {
  float tmp[3];    

  for (int i = 0; i < 3; i++) {
    std::copy(p_v, p_v + 3, tmp);
    p_v[(i+1)%3] = rc[i]*tmp[(i+1)%3] - rs[i]*tmp[(i+2)%3];
    p_v[(i+2)%3] = rc[i]*tmp[(i+2)%3] + rs[i]*tmp[(i+1)%3];
  }  
}

void tledMovingRigidContactSurfaceCPU::ApplyRotationsAndTranslationToNodeList(float *p_nodes, const float rc[], const float rs[], const float tPCor[]) {
  using namespace tledVectorArithmetic;

  float const *pc_n0 = this->GetAllNodeCoordinates0();

  for (float *p_n = p_nodes; p_n < p_nodes + 3*this->GetNumberOfNodes(); p_n += 3, pc_n0 += 3) {
    Sub(p_n, pc_n0, this->GetRotationCentre());
    this->RotateVector(p_n, rc, rs);
    Add(p_n, p_n, tPCor);
  }
}

void tledMovingRigidContactSurfaceCPU::SetAllNodeCoordinates0(const float nodes[]) {
  if ((int)m_Nodes0.size() != 3*this->GetNumberOfNodes()) m_Nodes0.resize(3*this->GetNumberOfNodes());
  std::copy(nodes, nodes + 3*this->GetNumberOfNodes(), m_Nodes0.begin());
}

void tledMovingRigidContactSurfaceCPU::SetAllOldNodeCoordinates(const float nodes[]) {
  if ((int)m_OldNodeCoordinates.size() != 3*this->GetNumberOfNodes()) m_OldNodeCoordinates.resize(3*this->GetNumberOfNodes());
  std::copy(nodes, nodes + 3*this->GetNumberOfNodes(), m_OldNodeCoordinates.begin());
}

void tledMovingRigidContactSurfaceCPU::ApplyTranslationToNodeList(float *p_nodes, const float t[]) {
  float const *pc_n0 = this->GetAllNodeCoordinates0();

  for (float *p_n = p_nodes; p_n < p_nodes + 3*this->GetNumberOfNodes(); p_n += 3, pc_n0 += 3) tledVectorArithmetic::Add(p_n, pc_n0, t);
}

void tledMovingRigidContactSurfaceCPU::ApplyRotationsToNormals(float *p_normals, const float n0s[], const int numNorms, const float rc[], const float rs[]) {
  std::copy(n0s, n0s + 3*numNorms, p_normals);
  for (float *p_n = p_normals; p_n < p_normals + 3*numNorms; p_n += 3) this->RotateVector(p_n, rc, rs);
}

void tledMovingRigidContactSurfaceCPU::SetNumberOfNodes(const int numNodes) {
  tledRigidContactSurfaceCPU::SetNumberOfNodes(numNodes);
  this->m_NodeNormals0.resize(3*numNodes);
}

void tledMovingRigidContactSurfaceCPU::SetNumberOfFacets(const int numFacets) {
  tledRigidContactSurfaceCPU::SetNumberOfFacets(numFacets);
  m_FacetNormals0.resize(3*numFacets);
  m_OldFacetNormals.resize(3*numFacets);
}

tledMovingRigidContactSurfaceCPU* tledMovingRigidContactSurfaceCPU::CreateSurface(const std::string &type) {
  tledMovingRigidContactSurfaceCPU *p_surf = NULL;

  if (type == "T3") {
    p_surf = new tledMovingRigidContactSurfaceT3CPU();
  } else if (type == "Q4") {
    p_surf = new tledMovingRigidContactSurfaceQ4CPU();
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Surface type \"" << type << "\" not recognised.");
  }

  return p_surf;
} 
