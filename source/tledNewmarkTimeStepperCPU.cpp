// =========================================================================
// File:       tledNewmarkTimeStepperCPU.cpp
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

#include "tledNewmarkTimeStepperCPU.h"

#include <algorithm>
#include <cassert>

tledNewmarkTimeStepperCPU::tledNewmarkTimeStepperCPU(const int numNodes, const float dt, const float alpha, const float M[]) : Superclass(numNodes, dt), mpc_M(M), mc_Alpha(alpha) {
  mp_PreviousAccelerations = new float[3*numNodes];
  mp_PreviousVelocities = new float[3*numNodes];

  std::fill(mp_PreviousAccelerations, mp_PreviousAccelerations + 3*numNodes, 0.0f);
  std::fill(mp_PreviousVelocities, mp_PreviousVelocities + 3*numNodes, 0.0f);
}

tledNewmarkTimeStepperCPU::~tledNewmarkTimeStepperCPU() {
  delete[] mp_PreviousAccelerations;
  delete[] mp_PreviousVelocities;
}

void tledNewmarkTimeStepperCPU::EvolveDisplacements(const float effectiveF[]) {
  for (int i = 0; i < 3*this->GetNumberOfNodes(); i++) {
    this->GetCurrentAccelerations()[i] = (effectiveF[i]/mpc_M[i] - mc_Alpha*(this->GetPreviousVelocities()[i] + this->GetTimeStep()/2*mp_PreviousAccelerations[i]))/(1 + mc_Alpha*this->GetTimeStep()/2);
    this->GetCurrentVelocities()[i] = this->GetPreviousVelocities()[i] + this->GetTimeStep()/2*(this->GetCurrentAccelerations()[i] + this->GetPreviousAccelerations()[i]);
    this->GetNextDisplacements()[i] = this->GetCurrentDisplacements()[i] + this->GetTimeStep()*(this->GetCurrentVelocities()[i] + this->GetTimeStep()/2*this->GetCurrentAccelerations()[i]);
  }
}

void tledNewmarkTimeStepperCPU::FinishTimeStep() {
  Superclass::FinishTimeStep();
  mp_PreviousVelocities = this->ReplaceCurrentVelocityBuffer(mp_PreviousVelocities);
  assert(&this->GetPreviousVelocities()[0] != &this->GetCurrentVelocities()[0]);
  mp_PreviousAccelerations = this->ReplaceCurrentAccelerationBuffer(mp_PreviousAccelerations);
  assert(&this->GetPreviousAccelerations()[0] != &this->GetCurrentAccelerations()[0]);
}
