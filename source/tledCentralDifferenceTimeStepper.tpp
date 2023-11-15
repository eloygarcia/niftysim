// =========================================================================
// File:       tledCentralDifferenceTimeStepper.tpp
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
tledCentralDifferenceTimeStepper<TTimeStepperBase>::tledCentralDifferenceTimeStepper(const int numNodes) : Superclass(numNodes) {
  mp_PreviousDisplacements = new float[3*numNodes];
}

template <class TTimeStepperBase>
tledCentralDifferenceTimeStepper<TTimeStepperBase>::~tledCentralDifferenceTimeStepper() {
  delete[] mp_PreviousDisplacements;
}

template <class TTimeStepperBase>
float* tledCentralDifferenceTimeStepper<TTimeStepperBase>::GetCurrentDeltaDisplacements(float *p_DU) const {
  float *p_du = p_DU;
  float const *pc_uc = this->GetCurrentDisplacements(), *pc_up = this->GetPreviousDisplacements();

  while (p_du < p_DU + 3*this->GetNumberOfNodes()) *(p_du++) = *(pc_uc++) - *(pc_up++);

  return p_DU;
}

template <class TTimeStepperBase>
float* tledCentralDifferenceTimeStepper<TTimeStepperBase>::RotateDisplacementBuffers(float *p_inBuffer) {
  assert(p_inBuffer == NULL);
  mp_PreviousDisplacements = Superclass::RotateDisplacementBuffers(mp_PreviousDisplacements);
  assert(mp_PreviousDisplacements != NULL);

  return NULL;
}
