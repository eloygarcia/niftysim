// =========================================================================
// File:       tledVectorArithmetic.cpp
// Purpose:    Standard vector functions
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledVectorArithmetic.h"

#include <cassert>
#include <limits>
#include <algorithm>

/*
 * Solves the 2x2 sub equation after the first elimination.
 */
static inline void _Solve22(float &r_x, float &r_y, const float m22[][2], const float r2[2]) {
  float pm;

  if (std::fabs(m22[0][0]) > std::fabs(m22[1][0])) {
    pm = m22[1][0]/m22[0][0];
    r_y = (r2[1] - r2[0]*pm)/(m22[1][1] - m22[0][1]*pm);
    r_x = (r2[0] - m22[0][1]*r_y)/m22[0][0]; 
  } else {
    pm = m22[0][0]/m22[1][0];
    r_y = (r2[0] - r2[1]*pm)/(m22[0][1] - m22[1][1]*pm);
    r_x = (r2[1] - m22[1][1]*r_y)/m22[1][0];
  }
}

/*
 * PLU direct solver
 */
namespace tledVectorArithmetic {
  void SolveEquation3x3(float &r_x, float &r_y, float &r_z, const float a[], const float b[], const float c[], const float r[]) {
    float m22[2][2], r2[2], pm;
  
    if (std::fabs(a[0]) > std::fabs(a[1]) && std::fabs(a[0]) > std::fabs(a[2])) {
      pm = a[1]/a[0];
      m22[0][0] = b[1] - b[0]*pm;
      m22[0][1] = c[1] - c[0]*pm;
      r2[0] = r[1] - r[0]*pm;

      pm = a[2]/a[0];
      m22[1][0] = b[2] - b[0]*pm;
      m22[1][1] = c[2] - c[0]*pm;
      r2[1] = r[2] - r[0]*pm;	  

      _Solve22(r_y, r_z, m22, r2);
      r_x = (r[0] - b[0]*r_y - c[0]*r_z)/a[0];
    } else if (std::fabs(a[1]) > std::fabs(a[2])) {
      pm = a[0]/a[1];
      m22[0][0] = b[0] - b[1]*pm;
      m22[0][1] = c[0] - c[1]*pm;
      r2[0] = r[0] - r[1]*pm;

      pm = a[2]/a[1];
      m22[1][0] = b[2] - b[1]*pm;
      m22[1][1] = c[2] - c[1]*pm;
      r2[1] = r[2] - r[1]*pm;	      
      _Solve22(r_y, r_z, m22, r2);
      r_x = (r[1] - b[1]*r_y - c[1]*r_z)/a[1];
    } else {
      /*
       * Pivot: line 2
       */
      pm = a[0]/a[2];
      m22[0][0] = b[0] - b[2]*pm;
      m22[0][1] = c[0] - c[2]*pm;
      r2[0] = r[0] - r[2]*pm;

      pm = a[1]/a[2];
      m22[1][0] = b[1] - b[2]*pm;
      m22[1][1] = c[1] - c[2]*pm;
      r2[1] = r[1] - r[2]*pm;	      
      _Solve22(r_y, r_z, m22, r2);
      r_x = (r[2] - b[2]*r_y - c[2]*r_z)/a[2];
    }
  } 

  static float _FastAcosSqrt(const float x, const int sign) {
    static const float ys[] = {1.5708e+00f,
			       1.4120e+00f,
			       1.3453e+00f,
			       1.2934e+00f,
			       1.2490e+00f,
			       1.2094e+00f,
			       1.1731e+00f,
			       1.1392e+00f,
			       1.1071e+00f,
			       1.0766e+00f,
			       1.0472e+00f,
			       1.0188e+00f,
			       9.9116e-01f,
			       9.6418e-01f,
			       9.3774e-01f,
			       9.1174e-01f,
			       8.8608e-01f,
			       8.6068e-01f,
			       8.3548e-01f,
			       8.1041e-01f,
			       7.8540e-01f,
			       7.6039e-01f,
			       7.3531e-01f,
			       7.1011e-01f,
			       6.8472e-01f,
			       6.5906e-01f,
			       6.3305e-01f,
			       6.0661e-01f,
			       5.7964e-01f,
			       5.5202e-01f,
			       5.2360e-01f,
			       4.9422e-01f,
			       4.6365e-01f,
			       4.3161e-01f,
			       3.9770e-01f,
			       3.6137e-01f,
			       3.2175e-01f,
			       2.7741e-01f,
			       2.2551e-01f,
			       1.5878e-01f,
			       0.0f};
    static const int numSteps = sizeof(ys)/sizeof(float) - 1;
    static const float dx = 1.0f/numSteps;
    const float xNorm = std::min(std::max(x, 0.0f), 1.0f)/dx;
    const int baseInd = (int)xNorm;
    const float lx = xNorm - baseInd;

    float y;
    
    assert(xNorm - baseInd <= 1);
    assert(xNorm - baseInd >= 0);
    assert(baseInd <= numSteps);

    y = ys[baseInd]*(1 - lx) + ys[baseInd+1]*lx;
    if (sign < 0) y = tledPi - y;

    return y;
  }

  float ComputeAngleNormalised(const float a[], const float b[]) {
    using namespace tledVectorArithmetic;

    assert(std::isnan(Norm(a)) || std::isnan(Norm(b)) || (std::fabs(1 - Norm(a)) < 1e-3f && std::fabs(1 - Norm(b)) < 1e-3f));

    return std::acos(std::max(std::min(Dot(a, b), 1.0f), -1.0f));
  }

  float ComputeAngleFast(const float a[], const float b[]) {
    using namespace tledVectorArithmetic;

    float na, nb;

    na = Dot(a, a);
    nb = Dot(b, b);
    if (na == 0 || nb == 0) return 0;
    else {
      float dab;
      int sign;

      dab = Dot(a, b);
      if (std::isnan(dab)) return std::numeric_limits<float>::quiet_NaN();
      sign = (dab < 0? -1 : 1);
      
      return _FastAcosSqrt((dab*dab)/(na*nb), sign);
    }
  }

  float ComputeAngle(const float a[], const float b[]) {
    using namespace tledVectorArithmetic;

    float na, nb;

    na = Dot(a, a);
    nb = Dot(b, b);
    if (tledHelper::IsNumericallyZero(std::fabs(na), std::fabs(nb)) || tledHelper::IsNumericallyZero(std::fabs(nb), std::fabs(na))) return 0;
    else {
      float dab = Dot(a, b);
      
      return std::acos((1 - 2*(dab < 0.f))*std::sqrt(std::max(std::min(dab*dab/(na*nb), 1.0f), 0.f)));
    }
  }

  float* Interpolate(float *p_dst, const float a[], const float b[], const float t) {
    float tmp[3];

    assert(p_dst != a && p_dst != b);
    Add(p_dst, ScalarMul(tmp, a, (1 - t)), ScalarMul(p_dst, b, t));
  
    return p_dst;
  }

  /* From D. Eberly: 3D game engine design */
  float* ComputeQuaternionFromRotation(float *p_dst, const float R[]) {
    p_dst[3] = 0.f;
    for (int r = 0; r < 3; r++) p_dst[3] += R[3*r+r];

#ifndef NDEBUG
    for (int r = 0; r < 3; r++) {
      assert(R[3*r+r] + 1 >= 1e-3f);
    }
#endif

    if (p_dst[3] > 0) {
      p_dst[3] = std::sqrt(p_dst[3] + 1)/2;

      p_dst[0] = (R[1*3+2] - R[2*3+1])/(4*p_dst[3]);
      p_dst[1] = (R[2*3+0] - R[0*3+2])/(4*p_dst[3]);
      p_dst[2] = (R[0*3+1] - R[1*3+0])/(4*p_dst[3]);
    } else {
      if (R[0*3+0] > R[1*3+1] && R[0*3+0] > R[2*3+2]) {
    	p_dst[0] = std::sqrt(2*R[0] - p_dst[3] + 1)/2;
    	p_dst[1] = (R[0*3+1] + R[1*3+0])/(4*p_dst[0]);
    	p_dst[2] = (R[0*3+2] + R[2*3+0])/(4*p_dst[0]);
    	p_dst[3] = (R[1*3+2] - R[2*3+1])/(4*p_dst[0]);
      } else if (R[1*3+1] > R[2*3+2]) {
    	p_dst[1] = std::sqrt(2*R[1*3+1] - p_dst[3] + 1)/2;
    	p_dst[0] = (R[0*3+1] + R[1*3+0])/(4*p_dst[1]);
    	p_dst[2] = (R[1*3+2] + R[2*3+1])/(4*p_dst[1]);
    	p_dst[3] = (R[2*3+0] - R[0*3+2])/(4*p_dst[1]);
      } else {
    	assert(R[2*3+2] >= R[1*3+1] && R[2*3+2] >= R[0*3+0]);
    	p_dst[2] = std::sqrt(2*R[2*3+2] - p_dst[3] + 1)/2;
	
    	p_dst[0] = (R[0*3+2] + R[2*3+0])/(4*p_dst[2]);
    	p_dst[1] = (R[1*3+2] + R[2*3+1])/(4*p_dst[2]);
    	p_dst[3] = (R[0*3+1] - R[1*3+0])/(4*p_dst[2]);
      }
    }

    return p_dst;
  }

  /* From D. Eberly: 3D game engine design */
  float* ComputeRotationFromQuaternion(float *p_dst, const float q[]) {
    float sqrs[3][4];

    for (int r = 0; r < 3; r++) for (int c = r; c < 4; c++) {
	sqrs[r][c] = q[r]*q[c];
      }

    p_dst[0*3+0] = 1.f - 2*(sqrs[1][1] + sqrs[2][2]), p_dst[1*3+0] = 2*(sqrs[0][1] - sqrs[2][3]), p_dst[2*3+0] = 2*(sqrs[0][2] + sqrs[1][3]);
    p_dst[0*3+1] = 2*(sqrs[0][1] + sqrs[2][3]), p_dst[1*3+1] = 1.f - 2*(sqrs[0][0] + sqrs[2][2]), p_dst[2*3+1] = 2*(sqrs[1][2] - sqrs[0][3]);
    p_dst[0*3+2] = 2*(sqrs[0][2] - sqrs[1][3]), p_dst[1*3+2] = 2*(sqrs[1][2] + sqrs[0][3]), p_dst[2*3+2] = 1.f - 2*(sqrs[0][0] + sqrs[1][1]);

#ifndef NDEBUG
    for (int c = 0; c < 3; c++) {
      float m = 0.f;

      for (int r = 0; r < 3; r++) m += p_dst[r*3+c]*p_dst[r*3+c];
      assert(std::fabs(1 - std::fabs(m)) < 1e-3f);
    }
#endif

    return p_dst;
  }

  float* OrthogonaliseNormalised(float *p_dst, const float a[], const float b[], const bool doNormalise) {
    float aProj[3], obNorm;

    ScalarMul(aProj, a, Dot(a, b));
    Sub(p_dst, b, aProj);

    if (doNormalise) {
      if ((obNorm = Norm(p_dst)) > 1e-6f) {
	ScalarDiv(p_dst, obNorm);
      }
    }

    return p_dst;
  }
}

