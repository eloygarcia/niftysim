// =========================================================================
// File:       tledMembraneMaterialLinear.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include <cstdlib>
#include <iostream>

#include "tledMembraneMaterialLinear.h"

void tledMembraneMaterialLinear::ComputeStress(float *p_stressOut, const float strainIn[]) const {
  p_stressOut[0] = strainIn[0] + m_Nu*strainIn[1];
  p_stressOut[1] = m_Nu*strainIn[0] + strainIn[1];
  p_stressOut[2] = (1 - m_Nu)/2*strainIn[2];

  for (int cInd = 0; cInd < 3; cInd++) p_stressOut[cInd] *= m_E/(1 - m_Nu*m_Nu);
}

void tledMembraneMaterialLinear::SetParameters(const std::vector<float> &elasticParameters) { 
  if (elasticParameters.size() != 2) {
    std::cerr << "Expected 2 parameters, received " << elasticParameters.size() << std::endl;
    std::abort();
  }

  m_E = elasticParameters[0], m_Nu = elasticParameters[1]; 
}
