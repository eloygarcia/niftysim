// =========================================================================
// File:       tledNewmarkTimeStepper.tpp
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

template <class TTimeStepperBase>
tledNewmarkTimeStepper<TTimeStepperBase>::tledNewmarkTimeStepper(const int numNodes, const float dt) : Superclass(numNodes), mc_Dt(dt) {
  mp_CurrentVelocities = new float[3*numNodes];
  mp_CurrentAccelerations = new float[3*numNodes];
}

template <class TTimeStepperBase>
tledNewmarkTimeStepper<TTimeStepperBase>::~tledNewmarkTimeStepper() {
  delete[] mp_CurrentAccelerations;
  delete[] mp_CurrentVelocities;
}


template <class TTimeStepperBase>
float* tledNewmarkTimeStepper<TTimeStepperBase>::GetCurrentDeltaDisplacements(float *p_DU) const {
  float *p_du = p_DU;
  float const *pc_vc = this->GetCurrentVelocities(), *pc_ac = this->GetCurrentAccelerations();

  while (p_du < p_DU + 3*this->GetNumberOfNodes()) *(p_du++) = this->GetTimeStep()*(*(pc_vc++) + this->GetTimeStep()/2**(pc_ac++));

  return p_DU;
}

template <class TTimeStepperBase>
float* tledNewmarkTimeStepper<TTimeStepperBase>::ReplaceCurrentVelocityBuffer(float *p_newBuffer) {
  float *p_old = mp_CurrentVelocities;

  mp_CurrentVelocities = p_newBuffer;

  return p_old;
}

template <class TTimeStepperBase>
float* tledNewmarkTimeStepper<TTimeStepperBase>::ReplaceCurrentAccelerationBuffer(float *p_newBuffer) {
  float *p_old = mp_CurrentAccelerations;

  mp_CurrentAccelerations = p_newBuffer;

  return p_old;
}
