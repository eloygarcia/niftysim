// =========================================================================
// File:       tledContactSurfaceCPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBaseSurface>
float* tledContactSurfaceImplCPU<TBaseSurface>::ComputeShapeValues(float *p_buffer, const float xi, const float eta) {
  if (Facet::NumberOfVertices == 3) {
    p_buffer[0] = 1.0f - xi - eta;
    p_buffer[1] = xi;
    p_buffer[2] = eta;
  } else {
    p_buffer[0] = (1.0f - xi)*(1.0f - eta);
    p_buffer[1] = xi*(1.0f - eta);
    p_buffer[1] = xi*eta;
    p_buffer[1] = (1.0f - xi)*eta;
  }

  return p_buffer;
}

template <class TBaseSurface>
bool tledContactSurfaceImplCPU<TBaseSurface>::ProjectOntoFacet(float *p_xi, const float x[], const float projOp[]) {
  using namespace tledVectorArithmetic;

  if (Facet::NumberOfVertices == 3) {
    static const float facetMargin = 1e-3f;
    const float *R0 = projOp;
    const float *R1 = R0 + 4;
    const float *R2 = R0 + 8;

    float xi[3];

    /*
     * Collision detection based on "double *collide (double **R, double *j, double *k)" from 
     * http://gpwiki.org/index.php/Ray_Triangle_Collision
     * (segment-version of Moller-Trumbore ray/triangle algorithm)
     */
    xi[0] = Dot(R0, x) + R0[3];
    xi[1] = Dot(R1, x) + R1[3];
    if (xi[1] < -facetMargin || xi[0] < -facetMargin || xi[0] + xi[1] > 1 + facetMargin) return false;
    else {
      float sum;

      p_xi[0] = std::max(xi[0], 0.0f);
      p_xi[1] = std::max(xi[1], 0.0f);

      if ((sum = p_xi[0] + p_xi[1]) > 1) {
	p_xi[0] /= sum;
	p_xi[1] /= sum;
      }

      p_xi[2] = Dot(R2, x) + R2[3];

      return true;
    }
  } else {
    tledFatalNotYetImplementedError;

    return false;
  }
}
