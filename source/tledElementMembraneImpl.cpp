// =========================================================================
// File:       tledElementMembraneImpl.cpp
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

#include "tledElementMembraneImpl.h"
#include "tledVectorArithmetic.h"

template <>
void tledElementMembraneImpl<3>::ComputeElementBasis(float (*p_basis)[3], float &r_area, const float X[]) const {
  using namespace tledVectorArithmetic;

  for (int e = 0; e < 2; e++) Sub(p_basis[e], X + 3*m_FacetVertexIndices[e+1], X + 3*m_FacetVertexIndices[0]);
  r_area = Norm(Cross(p_basis[2], p_basis[0], p_basis[1]))/2;
  ScalarDiv(p_basis[0], Norm(p_basis[0]));
  ScalarDiv(p_basis[2], r_area*2);
  ScalarDiv(p_basis[1], Norm(Cross(p_basis[1], p_basis[2], p_basis[0])));
}
