// =========================================================================
// File:       tledShellMaterialLinearThickPlateDecorator.tpp
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

template <class TMembraneMaterial>
void tledShellMaterialLinearThickPlateDecorator<TMembraneMaterial>::ComputeShearStress(float *p_stressOut, const float shearStrain[]) const {
  const float shearD = 5.0f/6.0f*this->GetPlateShearG()*this->GetThickness();

  p_stressOut[0] = shearD*shearStrain[0];
  p_stressOut[1] = shearD*shearStrain[1];
}
