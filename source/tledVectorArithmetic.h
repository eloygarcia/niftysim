// =========================================================================
// File:       tledVectorArithmetic.cpp
// Purpose:    Standard vector functions (CPU)
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
#ifndef tledVectorArithmetic_H
#define tledVectorArithmetic_H

#include "tledMatrixFunctions.h"
#include "tledHelper.h"

/**
 * \brief Namespace for generic CPU arithmetic operations
 * \ingroup helper
 */
namespace tledVectorArithmetic {
  __tled_inline float* Cross(float *p_cross, const float x[], const float y[]) {
    p_cross[0] = x[1]*y[2] - x[2]*y[1];
    p_cross[1] = x[2]*y[0] - x[0]*y[2];
    p_cross[2] = x[0]*y[1] - x[1]*y[0];

    return p_cross;
  } 

  __tled_inline float* Sub(float *p_r, const float x[], const float y[]) {
    return ::MatSubtract<float>(x, y, 3, 1, p_r);
  } 

  __tled_inline float* Add(float *p_r, const float x[], const float y[]) {
    return ::MatAdd<float>(x, y, 3, 1, p_r);
  } 

  __tled_inline float* ScalarDiv(float *p_x, const float a) {
    for (float *p_comp = p_x; p_comp < p_x + 3; p_comp++) *p_comp /= a;

    return p_x;
  } 

  __tled_inline float* ScalarDiv(float *p_xa, const float x[], const float a) {
    for (int cInd = 0; cInd < 3; cInd++) p_xa[cInd] = x[cInd]/a;

    return p_xa;
  } 

  __tled_inline float* ScalarMul(float *p_x, const float a) {
    for (float *p_comp = p_x; p_comp < p_x + 3; p_comp++) *p_comp *= a;

    return p_x;
  } 

  __tled_inline float* ScalarMul(float *p_ax, const float x[], const float a) {
    for (int cInd = 0; cInd < 3; cInd++) p_ax[cInd] = x[cInd]*a;

    return p_ax;
  } 

  /**
   * \brief 3x3 matrix-vector product
   */
  __tled_inline float* MatMul(float *p_Ax, const float A[][3], const float x[]) {
    for (int rInd = 0; rInd < 3; rInd++) {
      p_Ax[rInd] = A[rInd][0]*x[0];
      for (int cInd = 1; cInd < 3; cInd++) p_Ax[rInd] += A[rInd][cInd]*x[cInd];
    }

    return p_Ax;    
  }

  /**
   * \brief Euclidean scalar product
   */
  __tled_inline float Dot(const float x[], const float y[]) {
    return ::Dot<float>(x, 3, y);
  } 

  /**
   * \brief L2 norm
   */
  __tled_inline float Norm(const float x[]) {
    return std::sqrt(Dot(x, x));
  } 

  /**
   * \brief Linearly interpolates between \"a\" and \"b\".
   *
   * The destination buffer cannot overlap with either \"a\" or \"b\".
   */
  float* Interpolate(float *p_dst, const float a[], const float b[], const float t);

  /**
   * Safely solves the equation:
   * \f$\left(a|b|c\right)\cdot(s_x, s_y, s_z)^T = r\f$
   * Where a, b, c, r are 3-element vectors.
   */
  void SolveEquation3x3(float &r_sX, float &r_sY, float &r_sZ, const float a[], const float b[], const float c[], const float r[]);

  /**
   * \brief Computes the angle (in radians) between two normalised vectors.
   */
  float ComputeAngleNormalised(const float a[], const float b[]);

  /**
   * \brief Removes the co-linear with "a" component from "b"
   * 
   * 
   * "a" must be a (approximately) normalised vector, "b" is assumed to be of a similar magnitude.<br />
   * Destination pointer p_dst may point to the same memory location as "b".<br />
   * Optional normalisation of the output available via boolean flag.
   */
  float* OrthogonaliseNormalised(float *p_dst, const float a[], const float b[], const bool doNormalise = false);

  /**
   * \brief Computes the angle between two arbitrary-length vectors.
   */
  float ComputeAngle(const float a[], const float b[]);

  /**
   * \brief Computes the angle between two arbitrary-length vectors, uses a faster linear-interpolation based acos (error < 0.1).
   */
  float ComputeAngleFast(const float a[], const float b[]);

  /**
   * \brief Converts a 3x3 rotation matrix to a unit quaternion.
   * 
   * Destionation buffer size must be >= 4 elements.
   */
  float* ComputeQuaternionFromRotation(float *p_dst, const float R[]);

  /**
   * \brief Computes a 3x3 rotation matrix from a unit quaternion.
   */
  float* ComputeRotationFromQuaternion(float *p_dst, const float q[]);
}
#endif /* tledVectorArithmetic_H */
