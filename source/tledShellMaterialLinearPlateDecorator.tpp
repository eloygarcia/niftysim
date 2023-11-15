// =========================================================================
// File:       tledShellMaterialLinearPlateDecorator.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TMembraneMaterial>
void tledShellMaterialLinearPlateDecorator<TMembraneMaterial>::SetParameters(const std::vector<float> &all_prms) {
  std::vector<float> prms = all_prms;

  m_Nu = prms.back();
  prms.pop_back();
  m_E = prms.back();
  prms.pop_back();

  Superclass::SetParameters(prms);
}

template <class TMembraneMaterial>
void tledShellMaterialLinearPlateDecorator<TMembraneMaterial>::ComputeBendingResponse(float *p_stressOut, const float curvs[]) const {
  const float a = this->GetBendingE()*this->GetThickness()*this->GetThickness()*this->GetThickness()/(12*(1 - GetBendingNu()*GetBendingNu()));

  p_stressOut[0] = a*(curvs[0] + curvs[1]*this->GetBendingNu());
  p_stressOut[1] = a*(curvs[0]*this->GetBendingNu() + curvs[1]);
  p_stressOut[2] = a*curvs[2]*(1 - this->GetBendingNu())/2;
}

template <class TMembraneMaterial>
void tledShellMaterialLinearPlateDecorator<TMembraneMaterial>::ComputeStress(float *p_stressOut, const float strainCurvatures[]) const {
  Superclass::ComputeStress(p_stressOut, strainCurvatures);
  this->ComputeBendingResponse(p_stressOut + TMembraneMaterial::NumberOfStressComponents, strainCurvatures + TMembraneMaterial::NumberOfStrainComponents);
}

#ifdef _GPU_
template <class TMembraneMaterial>
void tledShellMaterialLinearPlateDecorator<TMembraneMaterial>::InitHostGPURepresentation(tledShellMaterial::GPUMaterial &r_sdst) const {
  GPUMaterial &r_dst = static_cast<GPUMaterial&>(r_sdst);

  Superclass::InitHostGPURepresentation(r_dst);
  r_dst.BendingE = m_E;
  r_dst.BendingNu = m_Nu;
} 

#endif
