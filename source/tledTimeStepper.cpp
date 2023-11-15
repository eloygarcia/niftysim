// =========================================================================
// File:       tledTimeStepper.cpp
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

#include "tledTimeStepper.h"

#include <algorithm>
#include <cassert>

float* tledTimeStepper::RotateDisplacementBuffers(float *p_inBuffer) {
  if (p_inBuffer == NULL) {
    std::iter_swap(&mp_CurrentDisplacements, &mp_NextDisplacements);

    return NULL;
  } else {
    float *p_oldCurrent = mp_CurrentDisplacements;

    mp_CurrentDisplacements = mp_NextDisplacements;
    mp_NextDisplacements = p_inBuffer;
    assert(p_oldCurrent != mp_CurrentDisplacements);
    assert(p_oldCurrent != mp_NextDisplacements);

    return p_oldCurrent;
  }
}

tledTimeStepper::tledTimeStepper(const int numNodes) : mc_NumNodes(numNodes) {
  mp_CurrentDisplacements = new float[numNodes*3];
  mp_NextDisplacements = new float[numNodes*3];
}

tledTimeStepper::~tledTimeStepper() {
  delete[] mp_CurrentDisplacements;
  delete[] mp_NextDisplacements;
}
