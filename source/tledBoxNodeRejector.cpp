// =========================================================================
// File:       tledBoxNodeRejector.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledBoxNodeRejector.h"
#include "tledModel.h"

bool tledBoxNodeRejector::DoReject(const int nodeIndex) const {
  const float *x = this->GetNode(nodeIndex);
  
  for (int axis = 0; axis < 3; axis++) {
    if (!(x[axis] >= this->GetMinBound(axis) && x[axis] <= this->GetMaxBound(axis))) return true;
  }
  
  return false;
}

void tledBoxNodeRejector::InitialiseFromXMLSpec(const XMLNode rootNode) {
  std::vector<float> bounds = GetXMLTextAsVector<float>(rootNode);

  if (bounds.size() != 6) {
    tledFatalError("Box boundary restrictor requires 6 floating-point parameters: MIN_X MAX_X ... MAX_Z");
  }

  for (int c = 0; c < 3; c++) {
    this->SetMinBound(c, bounds[2*c]);
    this->SetMaxBound(c, bounds[2*c+1]);
  }
}

