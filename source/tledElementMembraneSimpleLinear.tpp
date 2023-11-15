// =========================================================================
// File:       tledElementMembraneSimpleLinear.tpp
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

#ifdef _GPU_
template <const int t_numFacetVertices>
void tledElementMembraneSimpleLinear<t_numFacetVertices>::InitGPU(tledElementMembrane::GPUElement &r_bcDst) {
  GPUElement &r_dst = static_cast<GPUElement&>(r_bcDst);

  Superclass::InitGPU(r_dst);

  for (int aInd = 0; aInd < 3; aInd++) r_dst.ElementBasis[aInd] = make_float3(m_ElementBasis[aInd][0], m_ElementBasis[aInd][1], m_ElementBasis[aInd][2]);
  for (int r = 0; r < 2; r++) {
    r_dst.PhiInvT[r].x = m_PhiInvT[r][0];
    r_dst.PhiInvT[r].y = m_PhiInvT[r][1];
  }
}
#endif
