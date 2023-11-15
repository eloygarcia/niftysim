// =========================================================================
// File:       tledTimeStepperCPU.cpp
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

#include "tledTimeStepperCPU.h"

#include <algorithm>

tledTimeStepperCPU::tledTimeStepperCPU(const int numNodes) : Superclass(numNodes) {
  std::fill(this->GetCurrentDisplacements(), this->GetCurrentDisplacements() + 3*numNodes, 0.0f);
}

void tledTimeStepperCPU::SetCurrentDisplacements(const float U[]) {
  std::copy(U, U + 3*this->GetNumberOfNodes(), this->GetCurrentDisplacements());
}

void tledTimeStepperCPU::FinishTimeStep() {
  this->RotateDisplacementBuffers(NULL);
}
