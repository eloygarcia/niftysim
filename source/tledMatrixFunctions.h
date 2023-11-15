// =========================================================================
// File:       tledMatrixFunctions.h
// Purpose:    Utility functions for matrices
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    March 2007
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledMatrixFunctions_H
#define tledMatrixFunctions_H

#include "tledHelper.h"

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <cmath>

// ***** DECLARATIONS *****

// Utilities
template <class T> void MatPrint(const T* A, const int r, const int c);	// Print A
template <class T> void MatPrint(const std::vector<T>* A, const int r, const int c);
template <class T> void Swap(T* a, T* b); // Swap values in a and b

// std::vector operations
template <class T> T Dot(const T* A, const int n, const T* B); // returns dot(A,B), A & B are n-vectors
template <class T> T Norm(const T* A, const int n); // returns norm(A), A is an n-vector
template <class T> T Norm(const std::vector<T> A);

// Matrix operations
template <class T> void MatInv33(const T A[3][3], T R[3][3]);	// Inverse of A 3x3
template <class T> void MatDet33(const T A[3][3], T* R);	// Determinant of A 3x3
template <class T> void MatDiag33(const T A[3][3], T R[3]);	// Extract diagonal elements of A

#ifdef __CUDACC__
template <typename T> 
__device__ __host__
#else
template <typename T>
#endif
T* MatSubtract(const T* A, const T* B, const int r, const int c, T* R);	// Subtract B from A

#ifdef __CUDACC__
template <typename T> 
__device__ __host__
#else
template <typename T>
#endif
T* MatAdd(const T* A, const T* B, const int r, const int c, T* R);	// Add B to A

template <class T> void PermSym(int i, int j, int k, T* eijk);	// Permutation symbol
template <class T> void MatTrace(const T* A, const int r, T* R);	// Trace of A rxr
template <class T> void GaussJ(T* A, const int rA, const int cA, T* b, const int cb); // Solve the linear systems a*x = b, also converts a to inv(a)
template <class T> void MatInverse(T* A, const int rA); // Compute the inverse of matrix a

/** 
 * \brief Computes the transpose of a matrix 
 * \param p_matrix matrix to be transposed, R/W argument
 */
template <typename T> T* MatTranspose(T *p_matrix, const int numRows, const int numCols);

/**
 * Computes the SVD of a 3x3 matrix
 */
void SVD3x3(float *U, float *S, float *V, const float matrix[]);

/**
 * Eigenvalue decomposition for 3x3 SPD matrices (Jacobi algorithm). 
 * \param p_EV 3x3-element linear array for eigenvectors
 * \param p_E 3-element array for eigenvalues
 */
void SymmetricEVD3x3(float *p_EV, float *p_E, const float A[]);

/**
 * 3x3 QR Decomposition
 */
void QRD3x3(float *p_Q, float *p_R, const float A[]);

/**
 * \brief MxN QR Decomposition.
 *
 * Computes the QR decomposition of a general matrix A. It returns the <i>transpose</i> Q, not Q itself.
 * \param p_q (out) the destination for the orthogonal matrix Q. Note: this routine returns the transpose of Q.
 * \param p_r (out) the destination for the upper triangular matrix R can be the same as A, thus saving one copy operation.
 */
template <typename TScalar>
void QRDMxN(TScalar *p_q, TScalar *p_r, const TScalar A[], const int numRows, const int numCols);

#ifdef __CUDACC__
template <typename T> 
__device__ __host__
#else
template <typename T>
#endif
T* MatInverse22(T* iA, const T A[]); // Compute the inverse of 2x2 matrix a

#ifdef __CUDACC__
template <typename T> 
__device__ __host__
#else
template <typename T>
#endif
T MatDet22(const T A[]); // Compute the determinant of 2x2 matrix a

// Matrix multiplication
template <class T> void MatMult33T33(T A[3][3], T B[3][3], T R[3][3]);	// Product A'*B, A 3x3, B 3x3
template <class T> void MatMult43T43(T A[4][3], T B[4][3], T R[3][3]);	// Product A'*B, A 4x3, B 4x3
template <class T> void MatMult83T83(T A[8][3], T B[8][3], T R[3][3]);	// Product A'*B, A 8x3, B 8x3
template <class T> void MatMult84T83(T A[8][4], T B[8][3], T R[4][3]);	// Product A'*B, A 8x4, B 8x3
template <class T> void MatMult4343T(T A[4][3], T B[4][3], T R[4][4]);	// Product A*B', A 4x3, B 4x3
template <class T> void MatMult8383T(T A[8][3], T B[8][3], T R[8][8]);	// Product A*B', A 8x3, B 8x3
template <class T> void MatMult8484T(T A[8][4], T B[8][4], T R[8][8]);	// Product A*B', A 8x4, B 8x4
template <class T> void MatMult4333T(T A[4][3], T B[3][3], T R[4][3]); // Product A*B', A 4x3, B 3x3
template <class T> void MatMult8333T(T A[8][3], T B[3][3], T R[8][3]); // Product A*B', A 8x3, B 3x3
template <class T> void MatMult8883 (T A[8][8], T B[8][3], T R[8][3]); // Product A*B,  A 8x8, B 8x3
template <class T> void MatMult8884 (T A[8][8], T B[8][4], T R[8][4]); // Product A*B,  A 8x8, B 8x4
template <class T> void MatMult612T61(T A[6][12], T B[6], T R[12]);	// Product A'*B, A 6x12, B 6x1
template <class T> void MatMult624T61(T A[6][24], T B[6], T R[24]);	// Product A'*B, A 6x24, B 6x1

#ifdef __CUDACC__
template <typename T> 
__device__ __host__
#else
template <typename T>
#endif
T* MatMultScalar(const T* A, const int rA, const int cA, const T b, T* R);	// Product A*b, A rxc, b 1x1

template <class T> void MatMult3443(T A[3][4], T B[4][3], T R[3][3]);	// Product A*B, A 3x4, B 4x3
template <class T> void MatMult3333(T A[3][3], T B[3][3], T R[3][3]);	// Product A*B, A 3x3, B 3x3
template <class T> void MatMult3334(T A[3][3], T B[3][4], T R[3][4]);	// Product A*B, A 3x3, B 3x4
template <class T> void MatMult3883(T A[3][8], T B[8][3], T R[3][3]);	// Product A*B, A 3x8, B 8x3
template <class T> void MatMult3338(T A[3][3], T B[3][8], T R[3][8]);	// Product A*B, A 3x3, B 3x8
template <class T> void MatMult83T84(T A[8][3], T B[8][4], T R[3][4]);	// Product A'*B, A 8x3, B 8x4
template <class T> void MatMult8334 (T A[8][3], T B[3][4], T R[8][4]); // Product A*B,  A 8x3, B 3x4
template <class T> void MatMult8443 (T A[8][4], T B[4][3], T R[8][3]); // Product A*B,  A 8x4, B 4x3
template <class T> void MatMult3343T(T A[3][3], T B[4][3], T R[3][4]);	// Product A*B', A 3x3, B 4x3
template <class T> void MatMult3383T(T A[3][3], T B[8][3], T R[3][8]);	// Product A*B', A 3x3, B 8x3

template <class T> T* MatMultAB(const T* A, const int rA, const int cA, const T* B, const int rB, const int cB, T* R); // R = A*B, must have cA=rB, R will be rA-by-cB
template <class T> void MatMultAtB(T* A, const int rA, const int cA, T* B, const int rB, const int cB, T* R); // R = A'*B, must have rA=rB, R will be cA-by-cB
template <class T> void MatMultdAB(T* A, const int rA, T* B, const int rB, const int cB, T* R); // R = A*B, where A is a diagonal matrix and only diagonal components are stored (i.e. a std::vector)

// Interpolating splines
template <class T> void HermiteSpline(std::vector<T>* x, std::vector<T>* y, std::vector<T>* xx, std::vector<T>* yy); // Interpolate data (x,y) at points xx, return in yy
template <class T> void HermiteSpline(std::vector<T>* x, std::vector<T>* y, std::vector<T>* C3, std::vector<T>* C2, std::vector<T>* C1); // Compute interpolation polynomial coeffs for data (x,y), return in C1-C3
template <class T> void HermiteSpline(std::vector<T>* x, std::vector<T>* y, std::vector<T>* C3, std::vector<T>* C2, std::vector<T>* C1, T* xx, T* yy); // Given data (x,y) and precomputed coeffs C1-C3, interpolate at point xx, return in yy


// ***** DEFINITIONS *****

// Utilities
template <class T> void MatPrint(const T* A, const int r, const int c)
{
	int i,j;
	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			std::cout << *A << " ";
			A++;
		}
		std::cout << std::endl;
	}
}

template <class T> void MatPrint(std::vector<T>* A, int r, int c)
{
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			std::cout << (*A)[i*c+j] << " ";
		}
		std::cout << std::endl;
	}
}

template <class T> void Swap(T* a, T* b)
{
   T tmp = *a;
   *a = *b;
   *b = tmp;
}

// Vector operations
template <class T> T Dot(const T* A, const int n, const T* B)
{
   T R = 0;
   T const *a = A;
   T const *b = B;
   for (int i = 0; i < n; i++)
      R += (*a++)*(*b++);
   return R;
}

template <class T> T Norm(const T* A, const int n)
{
   T R = 0;

   for (T const *pc_a = A; pc_a < A + n; pc_a++) R += (*pc_a)*(*pc_a);

   return std::sqrt(R);
}

template <class T> T Norm(const std::vector<T> A)
{
   T R = 0;
   for (int i = 0; i < (int)A.size(); i++)
      R += A[i]*A[i];
   return std::sqrt(R);
}

// Matrix operations
template <class T> void MatInv33(const T A[3][3], T R[3][3])
{
	T detA;
	MatDet33(A,&detA);
	R[0][0] = (A[1][1]*A[2][2] - A[1][2]*A[2][1])/detA;
	R[0][1] = (A[0][2]*A[2][1] - A[0][1]*A[2][2])/detA;
	R[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1])/detA;
	R[1][0] = (A[1][2]*A[2][0] - A[1][0]*A[2][2])/detA;
	R[1][1] = (A[0][0]*A[2][2] - A[0][2]*A[2][0])/detA;
	R[1][2] = (A[0][2]*A[1][0] - A[0][0]*A[1][2])/detA;
	R[2][0] = (A[1][0]*A[2][1] - A[1][1]*A[2][0])/detA;
	R[2][1] = (A[0][1]*A[2][0] - A[0][0]*A[2][1])/detA;
	R[2][2] = (A[0][0]*A[1][1] - A[0][1]*A[1][0])/detA;
}

template <class T> void MatDet33(const T A[3][3], T* R)
{
	*R =  A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[1][0]*A[0][1]*A[2][2]
		+ A[1][0]*A[0][2]*A[2][1] + A[2][0]*A[0][1]*A[1][2] - A[2][0]*A[0][2]*A[1][1];
}

template <class T> void MatDiag33(const T A[3][3], T R[3])
{
	for (int i = 0; i < 3; i++)
	{
		R[i] = A[i][i];
	}
}

#ifdef __CUDACC__
template <typename T> 
__device__ __host__
#else
template <typename T>
#endif
T* MatSubtract(const T* A, const T* B, const int r, const int c, T* R) {
  T const *pc_a, *pc_b;
  T *p_r;

  for (pc_a = A, pc_b = B, p_r = R; p_r < R + r*c; p_r++, pc_a++, pc_b++) *p_r = *pc_a - *pc_b;

  return R;
}

#ifdef __CUDACC__
template <typename T> 
__device__ __host__
#else
template <typename T>
#endif
T* MatAdd(const T* A, const T* B, const int r, const int c, T* R) {
  T const *pc_a, *pc_b;
  T *p_r;

  for (pc_a = A, pc_b = B, p_r = R; p_r < R + r*c; p_r++, pc_a++, pc_b++) *p_r = *pc_a + *pc_b;

  return R;
}

template <class T> void PermSym(int i, int j, int k, T* eijk)
{
	// Check inputs
	if ((i < 1) || (i > 3) || (j < 1) || (j > 3) || (k < 1) || (k > 3))
	{
	  tledLogErrorStream(tledHelper::FatalError() << "Invalid values passed to PermSym: " << i << " " << j << " " << k);
	}

	T X[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	X[0][i-1] = 1;
	X[1][j-1] = 1;
	X[2][k-1] = 1;
	MatDet33(X,eijk);
}

template <class T> void MatTrace(T* A, int r, T* R)
{
	*R = 0;
	for (int i = 0; i < r; i++)
	{
		*R += *A;
		A += (r+1);
	}
}

// Matrix multiplication
template <class T> void MatMult33T33(T A[3][3], T B[3][3], T R[3][3])
{
	int i,j,k;
	memset(R,0,sizeof(T)*9);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
			{
				R[i][j] += A[k][i]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult43T43(T A[4][3], T B[4][3], T R[3][3])
{
	int i,j,k;
	memset(R,0,sizeof(T)*9);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 4; k++)
			{
				R[i][j] += A[k][i]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult83T83(T A[8][3], T B[8][3], T R[3][3])
{
	int i,j,k;
	memset(R,0,sizeof(T)*9);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 8; k++)
			{
				R[i][j] += A[k][i]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult84T83(T A[8][4], T B[8][3], T R[4][3])
{
	int i,j,k;
	memset(R,0,sizeof(T)*12);
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 8; k++)
			{
				R[i][j] += A[k][i]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult4343T(T A[4][3], T B[4][3], T R[4][4])
{
	int i,j,k;
	memset(R,0,sizeof(T)*16);
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < 3; k++)
			{
				R[i][j] += A[i][k]*B[j][k];
			}
		}
	}
}

template <class T> void MatMult8383T(T A[8][3], T B[8][3], T R[8][8])
{
	int i,j,k;
	memset(R,0,sizeof(T)*64);
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			for (k = 0; k < 3; k++)
			{
				R[i][j] += A[i][k]*B[j][k];
			}
		}
	}
}

template <class T> void MatMult8484T(T A[8][4], T B[8][4], T R[8][8])
{
	int i,j,k;
	memset(R,0,sizeof(T)*64);
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			for (k = 0; k < 4; k++)
			{
				R[i][j] += A[i][k]*B[j][k];
			}
		}
	}
}

template <class T> void MatMult4333T(T A[4][3], T B[3][3], T R[4][3])
{
	int i,j,k;
	memset(R,0,sizeof(T)*12);
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
			{
				R[i][j] += A[i][k]*B[j][k];
			}
		}
	}
}

template <class T> void MatMult8333T(T A[8][3], T B[3][3], T R[8][3])
{
	int i,j,k;
	memset(R,0,sizeof(T)*24);
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
			{
				R[i][j] += A[i][k]*B[j][k];
			}
		}
	}
}

template <class T> void MatMult8883 (T A[8][8], T B[8][3], T R[8][3])
{
	int i,j,k;
	memset(R,0,sizeof(T)*24);
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 8; k++)
			{
				R[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult8884 (T A[8][8], T B[8][4], T R[8][4])
{
	int i,j,k;
	memset(R,0,sizeof(T)*32);
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < 8; k++)
			{
				R[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult612T61(T A[6][12], T B[6], T R[12])
{
	int i,j;
	memset(R,0,sizeof(T)*12);
	for (i = 0; i < 12; i++)
	{
		for (j = 0; j < 6; j++)
		{
			R[i] += A[j][i]*B[j];
		}
	}
}

template <class T> void MatMult624T61(T A[6][24], T B[6], T R[24])
{
	int i,j;
	memset(R,0,sizeof(T)*24);
	for (i = 0; i < 24; i++)
	{
		for (j = 0; j < 6; j++)
		{
			R[i] += A[j][i]*B[j];
		}
	}
}

#ifdef __CUDACC__
template <typename T> 
__device__ __host__
#else
template <typename T>
#endif
T* MatMultScalar(const T* A, const int rA, const int cA, const T b, T* R)
{
  T const *pc_a;
  T *p_r;

  for (pc_a = A, p_r = R; pc_a < A + rA*cA; pc_a++, p_r++) {
    *p_r = (*pc_a)*b;
  }

  return R;
}

template <class T> void MatMult3443 (T A[3][4], T B[4][3], T R[3][3])
{
	int i,j,k;
	memset(R,0,sizeof(T)*12);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 4; k++)
			{
				R[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult3333 (T A[3][3], T B[3][3], T R[3][3])
{
	int i,j,k;
	memset(R,0,sizeof(T)*9);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
			{
				R[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult3334 (T A[3][3], T B[3][4], T R[3][4])
{
	int i,j,k;
	memset(R,0,sizeof(T)*12);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < 3; k++)
			{
				R[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult3883 (T A[3][8], T B[8][3], T R[3][3])
{
	int i,j,k;
	memset(R,0,sizeof(T)*24);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 8; k++)
			{
				R[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult3338 (T A[3][3], T B[3][8], T R[3][8])
{
	int i,j,k;
	memset(R,0,sizeof(T)*24);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 8; j++)
		{
			for (k = 0; k < 3; k++)
			{
				R[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult83T84(T A[8][3], T B[8][4], T R[3][4])
{
	int i,j,k;
	memset(R,0,sizeof(T)*12);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < 8; k++)
			{
				R[i][j] += A[k][i]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult8334 (T A[8][3], T B[3][4], T R[8][4])
{
	int i,j,k;
	memset(R,0,sizeof(T)*32);
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < 3; k++)
			{
				R[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult8443 (T A[8][4], T B[4][3], T R[8][3])
{
	int i,j,k;
	memset(R,0,sizeof(T)*24);
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 4; k++)
			{
				R[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

template <class T> void MatMult3343T(T A[3][3], T B[4][3], T R[3][4])
{
	int i,j,k;
	memset(R,0,sizeof(T)*12);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < 3; k++)
			{
				R[i][j] += A[i][k]*B[j][k];
			}
		}
	}
}

template <class T> void MatMult3383T(T A[3][3], T B[8][3], T R[3][8])
{
	int i,j,k;
	memset(R,0,sizeof(T)*24);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 8; j++)
		{
			for (k = 0; k < 3; k++)
			{
				R[i][j] += A[i][k]*B[j][k];
			}
		}
	}
}

template <class T> T* MatMultAB(const T* A, const int rA, const int cA, const T* B, const int rB, const int cB, T* R)
{
   const T* pA = NULL; const T* pB = NULL; T* pR = NULL;
   memset(R,0,sizeof(T)*rA*cB);
   for (int i = 0; i < cB; i++) // Columns of R
   {
      pA = A;
      for (int j = 0; j < rA; j++) // Rows of R
      {
         pR = R + i + j*cB;
         pB = B + i;
         for (int k = 0; k < cA; k++)
         {
            *pR += (*pA)*(*pB);
            pA++;
            pB += cB;
         }
      }
   }

   return R;
}

template <class T> void MatMultAtB(T* A, const int rA, const int cA, T* B, const int rB, const int cB, T* R)
{
   T* pA = NULL; T* pB = NULL; T* pR = NULL;
   memset(R,0,sizeof(T)*cA*cB);
   for (int i = 0; i < cA; i++) // Rows of R
   {
      for (int j = 0; j < cB; j++) // Columns of R
      {
         pR = R + j + i*cB;
         pA = A + i;
         pB = B + j;
         for (int k = 0; k < rA; k++)
         {
            *pR += (*pA)*(*pB);
            pA += cA;
            pB += cB;
         }
      }
   }
}

template <class T> void MatMultdAB(T* A, const int rA, T* B, const int rB, const int cB, T* R)
{
   T* pA = NULL; T* pB = NULL; T* pR = NULL;
   memset(R,0,sizeof(T)*rA*cB);
   pA = A;
   pB = B;
   for (int i = 0; i < rA; i++)  // Rows of R
   {
      for (int j = 0; j < cB; j++) // Columns of R
      {
         pR = R + i*cB + j;
         *pR = *pA*(*pB);
         pB++;
      }
      pA++;
   }
}

// Interpolating splines
template <class T> void HermiteSpline(std::vector<T>* x, std::vector<T>* y, std::vector<T>* xx, std::vector<T>* yy)
{
   // Check inputs
   int n = (int)(x->size());
   if (n < 2) {
     tledLogErrorStream(tledHelper::Warning() << "Must provide 2 or more data points");
   }
   if (y->size() != n) {
     tledLogErrorStream(tledHelper::Warning() << "Inputs x and y must be same size");
   }
   int nn = (int)(xx->size());
   if ( ((*x)[0] > (*xx)[0]) || ((*x)[n-1] < (*xx)[nn-1]) ) {
     tledLogErrorStream(tledHelper::Warning() << "Interpolation points xx are outside the bounds of data points x: extrapolation not allowed");
   }
   for (int i = 1; i < (int)(x->size()); i++)	// Check x is monotonic
   {
     if ((*x)[i] < (*x)[i-1]) {
       tledLogErrorStream(tledHelper::Warning() << "Input x must be monotonic");
     }
   }

   // Compute endpoint tangents
   std::vector<T> m(n,0);
   m[0] = ((*y)[1] - (*y)[0])/(2*((*x)[1] - (*x)[0]));
   m[n-1] = ((*y)[n-1] - (*y)[n-2])/(2*((*x)[n-1] - (*x)[n-2]));
   // Internal tangents
   if (n > 2)
   {
      for (int i = 1; i < n-1; i++)
      {
         m[i] = ((*y)[i+1] - (*y)[i])/(2*((*x)[i+1] - (*x)[i])) + ((*y)[i] - (*y)[i-1])/(2*((*x)[i] - (*x)[i-1]));
      }
   }

   // Find interval containing first point in xx
   int k = 0;
   bool found = false;
   while (!found)
   {
      if ((*x)[k+1] >= (*xx)[0])
         found = true;
      else
         k++;
   }
   T h = (*x)[k+1] - (*x)[k];
   T C3 = h*(m[k]+m[k+1]) + 2*((*y)[k]-(*y)[k+1]);
   T C2 = -h*(2*m[k]+m[k+1]) - 3*((*y)[k]-(*y)[k+1]);
   T C1 = h*m[k];
   T C0 = (*y)[k];

   // Loop over xx vals and compute corresponding yy
   for (int i = 0; i < nn; i++)
   {
      T t = ((*xx)[i] - (*x)[k])/h;
      (*yy)[i] = C3*t*t*t + C2*t*t + C1*t + C0;

      if (i < nn-1)
      {
         if ((*xx)[i+1] > (*x)[k+1])
         {
            k++;
            h = (*x)[k+1] - (*x)[k];
            C3 = h*(m[k]+m[k+1]) + 2*((*y)[k]-(*y)[k+1]);
            C2 = -h*(2*m[k]+m[k+1]) - 3*((*y)[k]-(*y)[k+1]);
            C1 = h*m[k];
            C0 = (*y)[k];
         }
      }
   }
}

template <class T> void HermiteSpline(std::vector<T>* x, std::vector<T>* y, std::vector<T>* C3, std::vector<T>* C2, std::vector<T>* C1)
{
   // Check inputs
   int n = (int)(x->size());
   if (n < 2) {
     tledLogErrorStream(tledHelper::Warning() << "Must provide 2 or more data points");
   }
   if (y->size() != n) {
     tledLogErrorStream(tledHelper::Warning() << "Inputs x and y must be same size");
   }
   for (int i = 1; i < (int)(x->size()); i++)	// Check x is monotonic
   {
     if ((*x)[i] < (*x)[i-1]) {
	tledLogErrorStream(tledHelper::Warning() << "Input x must be monotonic");
     }
   }

   // Compute endpoint tangents
   std::vector<T> m(n,0);
   m[0] = ((*y)[1] - (*y)[0])/(2*((*x)[1] - (*x)[0]));
   m[n-1] = ((*y)[n-1] - (*y)[n-2])/(2*((*x)[n-1] - (*x)[n-2]));
   // Internal tangents
   if (n > 2)
   {
      for (int i = 1; i < n-1; i++)
      {
         m[i] = ((*y)[i+1] - (*y)[i])/(2*((*x)[i+1] - (*x)[i])) + ((*y)[i] - (*y)[i-1])/(2*((*x)[i] - (*x)[i-1]));
      }
   }
   // Compute coeffs
   for (int k = 0; k < n-1; k++)
   {
      T h = (*x)[k+1] - (*x)[k];
      (*C3)[k] = h*(m[k]+m[k+1]) + 2*((*y)[k]-(*y)[k+1]);
      (*C2)[k] = -h*(2*m[k]+m[k+1]) - 3*((*y)[k]-(*y)[k+1]);
      (*C1)[k] = h*m[k];
   }

   // Clean up
   m.resize(0);
}

template <class T> void HermiteSpline(std::vector<T>* x, std::vector<T>* y, std::vector<T>* C3, std::vector<T>* C2, std::vector<T>* C1, T* xx, T* yy)
{
   // Check inputs
   int n = (int)(x->size());
   if (n < 2) {
     tledLogErrorStream(tledHelper::Warning() << "Must provide 2 or more data points");
   }
   if (y->size() != n) {
     tledLogErrorStream(tledHelper::Warning() << "Inputs x and y must be same size");
   }
   if ( ((*x)[0] > *xx) || ((*x)[n-1] < *xx) ) {
     tledLogErrorStream(tledHelper::Warning() << "Interpolation point xx is outside the bounds of data points x: extrapolation not allowed");
   }
   for (int i = 1; i < (int)(x->size()); i++)	// Check x is monotonic
   {
     if ((*x)[i] < (*x)[i-1]) {
       tledLogErrorStream(tledHelper::Warning() << "Input x must be monotonic");
     }
   }
   if ( (C3->size() != n-1) || (C2->size() != n-1) || (C1->size() != n-1) ) {
     tledLogErrorStream(tledHelper::Warning() << "Must provide n-1 interpolation coeffs");
   }

   // Find interval containing xx
   int k = 0;
   bool found = false;
   while (!found)
   {
      if ((*x)[k+1] >= *xx)
         found = true;
      else
         k++;
   }

   // Compute value
   T h = (*x)[k+1] - (*x)[k];
   T t = (*xx - (*x)[k])/h;
   *yy = (*C3)[k]*t*t*t + (*C2)[k]*t*t + (*C1)[k]*t + (*y)[k];
}

template <typename TScalar>
void QRDMxN(TScalar *p_q, TScalar *p_r, const TScalar a[], const int numRows, const int numCols) {
  std::vector<TScalar> tmpRRow(2*numCols);
  std::vector<TScalar> tmpQCol(2*numRows);

  std::fill(p_q, p_q + numRows*numRows, TScalar(0));
  for (int i = 0; i < numRows; i++) p_q[i*numRows+i] = TScalar(1);
  if (p_r != a) {
    std::copy(a, a + numRows*numCols, p_r);
  }

  /* Accumulating Q^T instead of Q -> row-major, convenience */
  for (int cInd = 0; cInd < numCols; cInd++) {
    for (int rInd = numRows - 1; rInd > cInd; rInd--) {
      TScalar killCos, killSin;

      assert(numRows > 0);

      {
	TScalar s = p_r[(rInd-1)*numCols+cInd], t = p_r[rInd*numCols+cInd];

	if (t + s == s) {
	  killCos = TScalar(1);
	  killSin = TScalar(0);
	} else {
	  if (std::fabs(t) > std::fabs(s)) {
	    TScalar theta = s/t;

	    killSin = TScalar(1)/std::sqrt(1 + theta*theta);
	    killCos = killSin*theta;
	  } else {
	    TScalar theta = t/s;

	    killCos = TScalar(1)/std::sqrt(1 + theta*theta);
	    killSin = killCos*theta;
	  }
	}
      }

      for (int multCInd = 0; multCInd < numCols; multCInd++) {
      	tmpRRow[multCInd] = p_r[(rInd-1)*numCols+multCInd]*killCos + p_r[rInd*numCols+multCInd]*killSin;
      }
      for (int multCInd = 0; multCInd < numCols; multCInd++) {
      	tmpRRow[numCols+multCInd] = -p_r[(rInd-1)*numCols+multCInd]*killSin + p_r[rInd*numCols+multCInd]*killCos;
      }
      std::copy(tmpRRow.begin(), tmpRRow.end(), p_r + (rInd - 1)*numCols);

      for (int multCInd = 0; multCInd < numRows; multCInd++) {
      	tmpQCol[multCInd] = p_q[(rInd-1)*numRows+multCInd]*killCos + p_q[rInd*numRows+multCInd]*killSin;
      }
      for (int multCInd = 0; multCInd < numRows; multCInd++) {
      	tmpQCol[numRows+multCInd] = -p_q[(rInd-1)*numRows+multCInd]*killSin + p_q[rInd*numRows+multCInd]*killCos;
      }
      std::copy(tmpQCol.begin(), tmpQCol.end(), p_q + (rInd - 1)*numRows);
    }
  }
}

template <typename TScalar>
void MatInverse(TScalar *p_A, const int numRowsCols) {
  std::vector<TScalar> qT(numRowsCols*numRowsCols), R(numRowsCols*numRowsCols), rInv(numRowsCols*numRowsCols, TScalar(0));

  QRDMxN(&qT.front(), &R.front(), p_A, numRowsCols, numRowsCols);
  for (int rInd = numRowsCols - 1; rInd >= 0; rInd--) {
    if (R[rInd*numRowsCols+rInd] == TScalar(0)) {
      tledNonFatalError("Matrix is singular aborting...");
      rInv[rInd*numRowsCols+rInd] = std::numeric_limits<TScalar>::quiet_NaN();
      return;
    } else rInv[rInd*numRowsCols+rInd] = TScalar(1)/R[rInd*numRowsCols+rInd];

    for (int cInd = rInd + 1; cInd < numRowsCols; cInd++) {
      TScalar acc = TScalar(0);

      for (int rInd2 = rInd; rInd2 < cInd; rInd2++) acc += rInv[rInd*numRowsCols+rInd2]*R[rInd2*numRowsCols+cInd];
      rInv[rInd*numRowsCols+cInd] = -acc/R[cInd*numRowsCols+cInd];
    }
  }

  MatMultAB(&rInv.front(), numRowsCols, numRowsCols, &qT.front(), numRowsCols, numRowsCols, p_A);
}

#ifdef __CUDACC__
template <typename T> 
__device__ __host__
#else
template <typename T>
#endif
T MatDet22(const T A[]) {
  return A[0*2+0]*A[1*2+1] - A[0*2+1]*A[1*2+0];
}

#ifdef __CUDACC__
template <typename T> 
__device__ __host__
#else
template <typename T>
#endif
T* MatInverse22(T *p_inv, const T A[]) {
  const T det = MatDet22(A);
  
  p_inv[0*2+0] = A[1*2+1]/det;
  p_inv[0*2+1] = -A[0*2+1]/det;
  p_inv[1*2+0] = -A[1*2+0]/det;
  p_inv[1*2+1] = A[0*2+0]/det;

  return p_inv;
}

template <typename T> 
T* MatTranspose(T *p_matrix, const int numRows, const int numCols) {
  for (int r = 0; r < numRows; r++) for (int c = r + 1; c < numCols; c++) std::iter_swap(p_matrix + r*3 + c, p_matrix + c*3 + r);

  return p_matrix;
}

#endif // tledMatrixFunctions_H
