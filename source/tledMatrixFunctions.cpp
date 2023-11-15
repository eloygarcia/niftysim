// =========================================================================
// File:       tledMatrixFunctions.cpp
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

#include "tledMatrixFunctions.h"
#include "tledVectorArithmetic.h"

#include <cmath>
#include <limits>
#include <cassert>
#include <algorithm>

static void _ApplyRotationLeft(float *p_M, const int i, const int j, const double c, const double s, const float M[]) {
  int remInd = -1;

  assert(i < j && i >= 0 && i < 3 && j >= 1 && j < 3);

  if (i == 0 && j == 1) remInd = 2;
  else if (i == 0 && j == 2) remInd = 1;
  else remInd = 0;

  for (int colInd = 0; colInd < 3; colInd++) {
    p_M[i*3+colInd] = (float)(c*M[3*i+colInd] + s*M[3*j+colInd]);
    p_M[j*3+colInd] = (float)(-s*M[3*i+colInd] + c*M[3*j+colInd]);
    p_M[remInd*3+colInd] = M[remInd*3+colInd];
  }  
}

static void _ApplyRotationRight(float *p_M, const int i, const int j, const double c, const double s, const float M[]) {
  int remInd = -1;

  assert(i < j && i >= 0 && i < 3 && j >= 1 && j < 3);

  if (i == 0 && j == 1) remInd = 2;
  else if (i == 0 && j == 2) remInd = 1;
  else remInd = 0;

  for (int rowInd = 0; rowInd < 3; rowInd++) {
    p_M[rowInd*3+i] = (float)(c*M[3*rowInd+i] + s*M[3*rowInd+j]);
    p_M[rowInd*3+j] = (float)(-s*M[3*rowInd+i] + c*M[3*rowInd+j]);
    p_M[rowInd*3+remInd] = M[rowInd*3+remInd];
  }  
}

void SymmetricEVD3x3(float *p_EV, float *p_E, const float A[]) {
  float S[9], SNew[9], EVNew[9];

  std::copy(A, A + 9, S);
  std::fill(p_EV, p_EV + 9, 0.0f);
  for (int k = 0; k < 3; k++) p_EV[3*k+k] = 1;

  for (int it = 0; it < 30; it++) {
    int maxRInd = -1, maxCInd = -1;
    float s, c, t;
    float pVal = -1, offDiag;

    offDiag = 0;
    for (int r = 0; r < 2; r++) for (int c = r + 1; c < 3; c++) {
	if (pVal < std::fabs(S[3*r+c])) {
	  pVal = std::fabs(S[3*r+c]);
	  maxRInd = r;
	  maxCInd = c;
	}
	offDiag += std::fabs(S[3*r+c]);
      }
    if (maxRInd < 0) return;
    
    if (std::fabs(S[0]) + std::fabs(S[4]) + std::fabs(S[8]) + offDiag == std::fabs(S[0]) + std::fabs(S[4]) + std::fabs(S[8])) break;

    {
      const float h = (S[maxRInd*3+maxRInd] - S[maxCInd*3+maxCInd]);
      const float theta  = 0.5f*h/S[maxRInd*3+maxCInd];

      t = 1/(std::fabs(theta) + std::sqrt(1 + theta*theta));
      if (theta < 0) t *= -1;

      c = 1/std::sqrt(1 + t*t);
      s = t*c;
    }

    _ApplyRotationLeft(SNew, maxRInd, maxCInd, c, s, S);    
    _ApplyRotationRight(S, maxRInd, maxCInd, c, s, SNew);

    _ApplyRotationRight(EVNew, maxRInd, maxCInd, c, s, p_EV);    
    std::copy(EVNew, EVNew + 9, p_EV);
  }

  for (int k = 0; k < 3; k++) p_E[k] = S[3*k+k];
}
  
static void _SwapSV(float *p_U, float *p_S, float *p_V, const int i, const int j) {
  std::iter_swap(p_S + i, p_S + j);
  for (int k = 0; k < 3; k++) {
    std::iter_swap(p_U + 3*k + i, p_U + 3*k + j);
    std::iter_swap(p_V + 3*k + i, p_V + 3*k + j);
  }
}

void SVD3x3(float *p_U, float *p_S, float *p_V, const float matrix[]) {
  float AA[3*3];

  for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) {
      AA[r*3+c] = 0;
      for (float const *pc_ar = matrix + c*3, *pc_al = matrix + r*3; pc_ar < matrix + c*3 + 3; pc_ar++, pc_al++) {
	  AA[r*3+c] += (*pc_al)*(*pc_ar);
	}
    }
  SymmetricEVD3x3(p_U, p_S, AA);

  for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) {
      AA[r*3+c] = 0;
      for (float const *pc_ar = matrix + c, *pc_al = matrix + r; pc_ar < matrix + c + 3*3; pc_ar += 3, pc_al += 3) {
	  AA[r*3+c] += (*pc_al)*(*pc_ar);
	}
    }
  SymmetricEVD3x3(p_V, p_S, AA);

  for (int d = 0; d < 3; d++) {
    if (p_S[d] < 0) p_S[d] = 0;
    else p_S[d] = std::sqrt(p_S[d]);
  }

  if (p_S[1] > p_S[0] && p_S[1] > p_S[2]) _SwapSV(p_U, p_S, p_V, 0, 1);
  else if (p_S[2] > p_S[0] && p_S[2] > p_S[1]) _SwapSV(p_U, p_S, p_V, 0, 2);
  if (p_S[2] > p_S[1]) _SwapSV(p_U, p_S, p_V, 1, 2);
}

void QRD3x3(float *p_Q, float *p_R, const float A[]) {
  using namespace tledVectorArithmetic;

  float tmp[3];

  for (int c = 0; c < 3; c++) p_Q[c] = A[3*c];
  p_R[0] = Norm(p_Q);
  p_R[3] = p_R[6] = 0;
  ScalarDiv(p_Q, p_R[0]);
  for (int c = 0; c < 3; c++) p_Q[3+c] = A[3*c+1];
  p_R[1] = Dot(p_Q + 3, p_Q);
  Sub(p_Q + 3, p_Q + 3, ScalarMul(tmp, p_Q, p_R[1]));
  p_R[4] = Norm(p_Q + 3);
  ScalarDiv(p_Q, p_R[4]);
  p_R[7] = 0;
  Cross(p_Q + 6, p_Q, p_Q + 3);
  for (int c = 0; c < 3; c++) tmp[c] = A[3*c+2];
  for (int c = 0; c < 3; c++) p_R[3*c+2] = Dot(tmp, p_Q + 3*c);
  MatTranspose(p_Q, 3, 3);
}
