// =========================================================================
// File:       tledSphereNodeRejector.cpp
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

#include "tledSphereNodeRejector.h"
#include "tledModel.h"
#include "tledVectorArithmetic.h"
#include "tledHelper.h"

bool tledSphereNodeRejector::DoReject(const int nodeIndex) const {
  using namespace tledVectorArithmetic;

  float d[3];
    
  return !(Norm(Sub(d, this->GetNode(nodeIndex), this->GetCentre())) < this->GetRadius());
}

void tledSphereNodeRejector::InitialiseFromXMLSpec(const XMLNode rootNode) {
  std::vector<float> spec = GetXMLTextAsVector<float>(rootNode);
  
  if (spec.size() != 4) {
    tledFatalError("Sphere boundary restrictor requires 4 floating-point parameters: RADIUS CENTRE_X CENTRE_Y CENTRE_Z");
  }

  this->SetRadius(spec[0]);
  this->SetCentre(&spec[1]);
}
