// =========================================================================
// File:       tledContactSolverCPU.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledContactSolverCPU.h"

void tledContactSolverCPU::ResetNormalResponse(ContactResponse &r_response) {
  std::fill(r_response.AccumulatedResponse, r_response.AccumulatedResponse + 3, 0.f);
  r_response.MaxProjection = -1.f;  
}

void tledContactSolverCPU::ResetFrictionResponse(ContactResponse &r_response) {
  std::fill(r_response.AccumulatedResponse, r_response.AccumulatedResponse + 3, 0.f);
  r_response.MaxProjection = -1.f;
}

void tledContactSolverCPU::SetContactForceBufferSize(const int numNodes, const bool doFriction) {
  m_NormalForces.resize(numNodes);  

  if (doFriction) {
    m_FrictionForces.resize(numNodes);  
    for (int i = 0; i < numNodes; i++) {
      this->ResetNormalResponse(i);
      this->ResetFrictionResponse(i);
    }
  } else {
    for (int i = 0; i < numNodes; i++) this->ResetNormalResponse(i);
  }
}
