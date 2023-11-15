// =========================================================================
// File:       tledElementShellBSTP1Edge.cu
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
#ifndef tledElementShellBSTP1Edge_CU
#define tledElementShellBSTP1Edge_CU

void tledElementShellBSTP1Edge::InitGPU(GPUElement &r_dst) {
#ifndef NDEBUG
  for (int i = 0; i < 2; i++) {
    r_dst.DeformationGradient[i] = make_float3(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
  }
  
  r_dst.NeighbourX0 = make_float3(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
  for (int n = 0; n < 4; n++) {
    r_dst.PatchShapeGradients[n] = make_float2(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
  }
#endif

  if (m_NumPatchNodes == 3) {
    r_dst.NeighbourNodeIndex = -1;   
    r_dst.EdgeNormalAndLength0.x = m_EdgeNormal[0];
    r_dst.EdgeNormalAndLength0.y = m_EdgeNormal[1];
    r_dst.EdgeNormalAndLength0.z = m_L0;
    r_dst.PhiN0 = make_float3(m_PhiN0[0], m_PhiN0[1], m_PhiN0[2]);
  } else {
    r_dst.NeighbourX0 = make_float3(*(mc_X0s + 3*m_EdgePatchNodeIndices[3]), *(mc_X0s + 3*m_EdgePatchNodeIndices[3] + 1), *(mc_X0s + 3*m_EdgePatchNodeIndices[3] + 2));
    r_dst.NeighbourPatchShapeGradient = make_float2(m_PatchShapeGradients[3][0], m_PatchShapeGradients[3][1]);
    r_dst.NeighbourNodeIndex = m_EdgePatchNodeIndices[3];
  }

  for (int n = 0; n < 3; n++) {
    r_dst.PatchShapeGradients[n].x = m_PatchShapeGradients[n][0];
    r_dst.PatchShapeGradients[n].y = m_PatchShapeGradients[n][1];
  }
}

#endif
