// =========================================================================
// File:       tledCentralDifferenceTimeStepperCPU.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledCentralDifferenceTimeStepperCPU.h"

#include <algorithm>

void tledCentralDifferenceTimeStepperCPU::FinishTimeStep() {
  this->RotateDisplacementBuffers(NULL);
}

void tledCentralDifferenceTimeStepperCPU::EvolveDisplacements(const float effectiveF[]) {
  for (int i = 0; i < 3*this->GetNumberOfNodes(); i++) {
    this->GetNextDisplacements()[i] = mp_A[i]*effectiveF[i] + mp_B[i]*this->GetCurrentDisplacements()[i] + mp_C[i]*this->GetPreviousDisplacements()[i];
  }
}

tledCentralDifferenceTimeStepperCPU::tledCentralDifferenceTimeStepperCPU(const int numNodes, const float dt, const float alpha, const float M[]) : Superclass(numNodes) {
  mp_A = new float[3*numNodes];
  mp_B = new float[3*numNodes];
  mp_C = new float[3*numNodes];

  std::fill(this->GetPreviousDisplacements(), this->GetPreviousDisplacements() + 3*numNodes, 0.0f);

  for (int i = 0; i < 3*numNodes; i++) {
    mp_A[i] = 1/(M[i]/dt/dt + alpha*M[i]/2/dt);
    mp_B[i] = 2*M[i]*mp_A[i]/dt/dt;
    mp_C[i] = alpha*M[i]*mp_A[i]/2/dt - mp_B[i]/2;
  }
}

tledCentralDifferenceTimeStepperCPU::~tledCentralDifferenceTimeStepperCPU() {
  delete[] mp_A;
  delete[] mp_B;
  delete[] mp_C;
}

