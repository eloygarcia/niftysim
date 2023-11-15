// =========================================================================
// File:       tledElementShellBSTP1.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   Cuda
// Created:    August 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledElementShellBSTP1_CU
#define tledElementShellBSTP1_CU

#include <limits>

#include "tledElementShellBSTP1.h"
#include "tledElementShellBSTP1Edge.cu"

void tledElementShellBSTP1::InitGPU(tledElementMembrane::GPUElement &r_bcDst) {
  GPUElement &r_dst = static_cast<GPUElement&>(r_bcDst);

  Superclass::InitGPU(r_dst);
#ifndef NDEBUG
  for (int i = 0; i < 2; i++) {
    r_dst.ContravariantBasis[i] = make_float3(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
  }
  r_dst.CurrentNormal = make_float3(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
  for (int i = 0; i < 3; i++) {
    r_dst.LinearShapeGradients[i] = make_float2(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
  }
#endif

  for (int e = 0; e < 3; e++) m_Edges[e].InitGPU(r_dst.EdgeElements[e]);
  for (int n = 0; n < 3; n++) {
    r_dst.X0s[n] = make_float3(mc_X0[3*this->GetFacetVertexIndices()[n]], mc_X0[3*this->GetFacetVertexIndices()[n]+1], mc_X0[3*this->GetFacetVertexIndices()[n]+2]);
    r_dst.LinearShapeGradients[n] = make_float2(GetLinearGradient(n)[0], GetLinearGradient(n)[1]);
  }
}

#endif
