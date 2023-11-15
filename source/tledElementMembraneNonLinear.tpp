// =========================================================================
// File:       tledElementMembraneNonLinear.tpp
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

#include <algorithm>
#include <limits>

template <const int t_numFacetVertices>
tledElementMembraneNonLinear<t_numFacetVertices>::tledElementMembraneNonLinear() {
  std::fill(m_LastStrain, m_LastStrain + 8, 0.f);  
#ifndef NDEBUG
  std::fill(&m_InitialJacobian[0][0], &m_InitialJacobian[0][0] + 6, std::numeric_limits<float>::quiet_NaN());
  std::fill(&m_InitialCauchyGreenTensor[0][0], &m_InitialCauchyGreenTensor[0][0] + 4, std::numeric_limits<float>::quiet_NaN());  
#endif
}

#ifdef _GPU_
template <const int t_numFacetVertices>
void tledElementMembraneNonLinear<t_numFacetVertices>::InitGPU(tledElementMembrane::GPUElement &r_bcDst) {
  GPUElement &r_dst = static_cast<GPUElement&>(r_bcDst);

  Superclass::InitGPU(r_dst);
  for (int i = 0; i < 2; i++) {
    r_dst.InitialJacobian[i] = make_float3(m_InitialJacobian[i][0], m_InitialJacobian[i][1], m_InitialJacobian[i][2]);
    r_dst.InitialCauchyGreenTensor[i] = make_float2(m_InitialCauchyGreenTensor[i][0], m_InitialCauchyGreenTensor[i][1]);
  }
}
#endif
