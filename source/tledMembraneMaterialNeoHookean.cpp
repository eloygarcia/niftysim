// =========================================================================
// File:       tledMembraneMaterialNeoHookean.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledMatrixFunctions.h"
#include "tledMembraneMaterialNeoHookean.h"

void tledMembraneMaterialNeoHookean::ComputeStress(float *p_stressOut, const float def[]) const {
  ComputeStress(p_stressOut, def, def + 2*2);
}

void tledMembraneMaterialNeoHookean::ComputeStress(float *p_stressOut, const float C0[], const float Cn[]) const {
  float psk[2*2], invCn[2*2];
  float vRatio;

  MatInverse22(psk, C0);
  MatInverse22(invCn, Cn);

  vRatio = MatDet22(C0)/MatDet22(Cn);
  for (int e = 0; e < 4; e++) psk[e] -= vRatio*invCn[e];

  MatMultScalar(psk, 2, 2, m_Mu*GetThickness(), p_stressOut);
}

void tledMembraneMaterialNeoHookean::SetParameters(const std::vector<float> &elasticParameters) {
  m_Mu = elasticParameters.front();
}
