// =========================================================================
// File:       tledNarrowConeSelfCollisionBVH.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TSurface, class TBV, class TAPI>
void tledNarrowConeSelfCollisionBVH<TSurface, TBV, TAPI>::ComputeSurfaceConePairwise(float &r_angle, float *p_axis, const float angle0, const float axis0[], const float angle1, const float axis1[]) {
  using namespace tledVectorArithmetic;

  float beta0, beta1, wc0, wc1;
  float axisNorm, beta, tmp[3], e;

  e = Dot(axis0, axis1);
  if (e > -0.99f) {
    beta = ComputeAngleNormalised(axis0, axis1);
    beta0 = (2*beta - (angle0 - angle1))/4;
    beta1 = beta - beta0;	  
    wc1 = (std::cos(beta1) - e*std::cos(beta0))/(1 - e*e);
    wc0 = std::cos(beta0) - wc1*e;
    assert(wc0 >= -1e-4f && wc1 >= -1e-4f);
    ScalarMul(p_axis, axis0, wc0);
    Add(p_axis, p_axis, ScalarMul(tmp, axis1, wc1));
    if ((axisNorm = Norm(p_axis)) < 1e-6) {
      /* Angle between Vol. axes ~Pi */
      r_angle = 2*BoundingVolume::GetVolinoThresholdAngle();
    } else {
      assert(axisNorm != axisNorm || axisNorm > 0);
      ScalarDiv(p_axis, axisNorm);
      assert(ComputeAngleNormalised(p_axis, axis0) < tledPi && ComputeAngleNormalised(p_axis, axis1) < tledPi);
      r_angle = 2.0001f*std::max(ComputeAngleNormalised(p_axis, axis0) + angle0/2, ComputeAngleNormalised(p_axis, axis1) + angle1/2);
    }  
  } else {
    r_angle = 2*BoundingVolume::GetVolinoThresholdAngle();
  }

  assert(r_angle >= 0);
}
