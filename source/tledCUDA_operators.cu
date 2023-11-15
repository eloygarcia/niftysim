// =========================================================================
// File:       tledCUDA_operators.cu
// Purpose:    Overloaded operators for CUDA vector types
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   cu
// Created:    May 2007
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledCUDA_operators_CU
#define tledCUDA_operators_CU

#ifdef _CUDA_5PLUS_SDK
#include <helper_math.h>
#endif

// Operator "+" ~~~~~~~~~~~~~~~~~~~~~~
#ifndef _CUDA_5PLUS_SDK
inline __device__ float3 operator+(const float3 a, const float3 b)
{
   return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
#endif

inline __device__ float3 operator+(const float3 &a, const float4 &b) {
  return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __device__ float2 operator+(const float2 &a, const float2 &b) {
  return make_float2(a.x + b.x, a.y + b.y);
}

#ifndef _CUDA_5PLUS_SDK
inline __device__ float3 operator+(const float3 a, const float b)
{
   return make_float3(a.x + b, a.y + b, a.z + b);
}

inline __device__ float3 operator+(const float a, const float3 b)
{
   return make_float3(a + b.x, a + b.y, a + b.z);
}
#endif

inline __device__ float3 operator+(const float3 a, const int b)
{
   return make_float3(a.x + (float)b, a.y + (float)b, a.z + (float)b);
}

inline __device__ float3 operator+(const int a, const float3 b)
{
   return make_float3((float)a + b.x, (float)a + b.y, (float)a + b.z);
}

#ifndef _CUDA_5PLUS_SDK
inline __device__ int3 operator+(const int3 a, const int b)
{
   return make_int3(a.x + b, a.y + b, a.z + b);
}

inline __device__ int3 operator+(const int a, const int3 b)
{
   return make_int3(a + b.x, a + b.y, a + b.z);
}

inline __device__ float4 operator+(const float4 a, const float4 b)
{
   return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

inline __device__ int4 operator+(const int4 a, const int b)
{
   return make_int4(a.x + b, a.y + b, a.z + b, a.w + b);
}
#endif

// Operator "+=" ~~~~~~~~~~~~~~~~~~~~~~
#ifndef _CUDA_5PLUS_SDK
inline __device__ void operator+=(float4 &a, const float4 &b) {
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  a.w += b.w;
}

inline __device__ void operator+=(float3 &a, const float3 &b) {
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
}
#endif

inline __device__ void operator+=(float4 &r_a, const float3 &b) {
  r_a.x += b.x;
  r_a.y += b.y;
  r_a.z += b.z;
}

inline __device__ void operator+=(float2 &r_a, const float2 &b) {
  r_a.x += b.x;
  r_a.y += b.y;
}

// Operator "-" ~~~~~~~~~~~~~~~~~~~~~~
#ifndef _CUDA_5PLUS_SDK
inline __device__ float3 operator-(const float3 a, const float3 b)
{
   return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __device__ float3 operator-(const float3 a, const float b)
{
   return make_float3(a.x - b, a.y - b, a.z - b);
}

inline __device__ float3 operator-(const float a, const float3 b)
{
   return make_float3(a - b.x, a - b.y, a - b.z);
}
#endif

inline __device__ float3 operator-(const float3 a, const int b)
{
   return make_float3(a.x - (float)b, a.y - (float)b, a.z - (float)b);
}

inline __device__ float3 operator-(const int a, const float3 b)
{
   return make_float3((float)a - b.x, (float)a - b.y, (float)a - b.z);
}

inline __device__ float3 operator-(const float3 a)
{
   return make_float3(-a.x, -a.y, -a.z);
}

#ifndef _CUDA_5PLUS_SDK
inline __device__ float4 operator-(const float4 a, const float4 b)
{
   return make_float4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}
#endif

inline __device__ float4 operator-(const float4 a)
{
   return make_float4(-a.x, -a.y, -a.z, -a.w);
}

// Operator "-=" ~~~~~~~~~~~~~~~~~~~~~~
#ifndef _CUDA_5PLUS_SDK
inline __device__ void operator-=(float4& r_a, const float4 &b) 
{
  r_a.x -= b.x;
  r_a.y -= b.y;
  r_a.z -= b.z;
  r_a.w -= b.w;
}

inline __device__ void operator-=(float3& r_a, const float3 &b)
{
  r_a.x -= b.x;
  r_a.y -= b.y;
  r_a.z -= b.z;
}
#endif

inline __device__ void operator-=(float4& r_a, const float3 &b) {
  r_a.x -= b.x;
  r_a.y -= b.y;
  r_a.z -= b.z;
}

// Operator "*=" ~~~~~~~~~~~~~~~~~~~~~~
#ifndef _CUDA_5PLUS_SDK
inline __device__ void operator*=(float3 &r_a, const float s) 
{
  r_a.x *= s;
  r_a.y *= s;
  r_a.z *= s;
}

inline __device__ void operator*=(float2 &r_a, const float s) 
{
  r_a.x *= s;
  r_a.y *= s;
}
#endif

// Operator "/=" ~~~~~~~~~~~~~~~~~~~~~~
#ifndef _CUDA_5PLUS_SDK
inline __device__ void operator/=(float3 &r_a, const float s) {
  r_a.x /= s;
  r_a.y /= s;
  r_a.z /= s;
}

inline __device__ void operator/=(float2 &r_a, const float s) {
  r_a.x /= s;
  r_a.y /= s;
}
#endif

// Operator "*" ~~~~~~~~~~~~~~~~~~~~~~
#ifndef _CUDA_5PLUS_SDK
inline __device__ float3 operator*(const float a, const float3 b)
{
   return make_float3(a*b.x, a*b.y, a*b.z);
}

inline __device__ float3 operator*(const float3 a, const float b)
{
   return make_float3(a.x*b, a.y*b, a.z*b);
}

inline __device__ float3 operator*(const float3 a, const float3 b)
{
   // Piece-wise multiplication (like .* in Matlab)
   return make_float3(a.x*b.x, a.y*b.y, a.z*b.z);
}

inline __device__ int3 operator*(const int3 a, const int3 b)
{
   // Piece-wise multiplication (like .* in Matlab)
   return make_int3(a.x*b.x, a.y*b.y, a.z*b.z);
}
#endif

inline __device__ float3 operator*(const float3 a, const int3 b)
{
   // Piece-wise multiplication (like .* in Matlab)
   return make_float3(a.x*((float)b.x), a.y*((float)b.y), a.z*((float)b.z));
}

inline __device__ float3 operator*(const int3 a, const float3 b)
{
   // Piece-wise multiplication (like .* in Matlab)
   return make_float3(((float)a.x)*b.x, ((float)a.y)*b.y, ((float)a.z)*b.z);
}

#ifndef _CUDA_5PLUS_SDK
inline __device__ float4 operator*(const float a, const float4 b)
{
   return make_float4(a*b.x, a*b.y, a*b.z, a*b.w);
}
#endif

inline __device__ float4 operator*(const float4 a, const int4 b)
{
   // Piece-wise multiplication (like .* in Matlab)
   return make_float4(a.x*((float)b.x), a.y*((float)b.y), a.z*((float)b.z), a.w*((float)b.w));
}

inline __device__ float4 operator*(const int4 a, const float4 b)
{
   // Piece-wise multiplication (like .* in Matlab)
   return make_float4(((float)a.x)*b.x, ((float)a.y)*b.y, ((float)a.z)*b.z, ((float)a.w)*b.w);
}

// Operator "/" ~~~~~~~~~~~~~~~~~~~~~~
#ifndef _CUDA_5PLUS_SDK
inline __device__ float4 operator/(const float4 a, const float b)
{
   // Scalar divide
   return make_float4((a.x)/b, (a.y)/b, (a.z)/b, (a.w)/b);
}

inline __device__ float3 operator/(const float3 a, const float b)
{
   // Scalar divide
   return make_float3((a.x)/b, (a.y)/b, (a.z)/b);
}
#endif

// Operator "%" ~~~~~~~~~~~~~~~~~~~~~~
inline __device__ int3 operator%(const int3 a, const int b)
{
   return make_int3(a.x%b, a.y%b, a.z%b);
}

inline __device__ int4 operator%(const int4 a, const int b)
{
   return make_int4(a.x%b, a.y%b, a.z%b, a.w%b);
}

// Operator "==" ~~~~~~~~~~~~~~~~~~~~~~
inline __device__ bool operator==(const int2 &a, const int2 &b) {
  return a.x == b.x && a.y == b.y;
}

inline __device__ bool operator!=(const int2 &a, const int2 &b) {
  return !(a == b);
}

// Matrix and vector functions
#ifndef _CUDA_5PLUS_SDK
inline __device__ float dot(const float3 a, const float3 b)
{
   return (a.x*b.x + a.y*b.y + a.z*b.z);
}

inline __device__ float dot(const float4 a, const float4 b)
{
   return (a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w);
}

inline __device__ float dot(const float2 &a, const float2 &b) {
  return a.x*b.x + a.y*b.y;
}
#endif

inline __device__ float dot(const float3 &a, const float4 &b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline __device__ float dot(const float4 &a, const float3 &b) {
  return dot(b, a);
}

inline __device__ float3 interpolate(const float3 &a, const float3 &b, const float t) {
  return a + t*(b - a);
}

#ifndef _CUDA_5PLUS_SDK
inline __device__ float3 cross(const float3 &a, const float3 &b) {
  return make_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
#endif

inline __device__ float norm(const float3 &a) {
  return sqrtf(dot(a, a));
}

inline __device__ float norm(const float2 &a) {
  return sqrtf(dot(a, a));
}

inline __device__ float sum(const float3 a)
{
   return a.x + a.y + a.z;
}

inline __device__ float sum(const float4 a)
{
   return a.x + a.y + a.z + a.w;
}

namespace tledCUDAMaths {
  /*
   * Solves the 2x2 sub equation after the first elimination.
   */
  static __device__ void _Solve2x2(float &r_x, float &r_y, const float m22[][2], const float r2[2]) {
    float pm;

    if (fabsf(m22[0][0]) > fabsf(m22[1][0])) {
      pm = m22[1][0]/m22[0][0];
      r_y = (r2[1] - r2[0]*pm)/(m22[1][1] - m22[0][1]*pm);
      r_x = (r2[0] - m22[0][1]*r_y)/m22[0][0]; 
    } else {
      pm = m22[0][0]/m22[1][0];
      r_y = (r2[0] - r2[1]*pm)/(m22[0][1] - m22[1][1]*pm);
      r_x = (r2[1] - m22[1][1]*r_y)/m22[1][0];
    }
  }

  inline __device__ void SolveEquation3x3(float &r_x, float &r_y, float &r_z, const float3 &a, const float3 &b, const float3 &c, const float3 &r) {
    float m22[2][2], r2[2], pm;
  
    if (fabsf(a.x) > fabsf(a.y) && fabsf(a.x) > fabsf(a.z)) {
      pm = a.y/a.x;
      m22[0][0] = b.y - b.x*pm;
      m22[0][1] = c.y - c.x*pm;
      r2[0] = r.y - r.x*pm;

      pm = a.z/a.x;
      m22[1][0] = b.z - b.x*pm;
      m22[1][1] = c.z - c.x*pm;
      r2[1] = r.z - r.x*pm;	  

      _Solve2x2(r_y, r_z, m22, r2);
      r_x = (r.x - b.x*r_y - c.x*r_z)/a.x;
    } else if (fabsf(a.y) > fabsf(a.z)) {
      pm = a.x/a.y;
      m22[0][0] = b.x - b.y*pm;
      m22[0][1] = c.x - c.y*pm;
      r2[0] = r.x - r.y*pm;

      pm = a.z/a.y;
      m22[1][0] = b.z - b.y*pm;
      m22[1][1] = c.z - c.y*pm;
      r2[1] = r.z - r.y*pm;	      
      _Solve2x2(r_y, r_z, m22, r2);
      r_x = (r.y - b.y*r_y - c.y*r_z)/a.y;
    } else {
      /*
       * Pivot: line 2
       */
      pm = a.x/a.z;
      m22[0][0] = b.x - b.z*pm;
      m22[0][1] = c.x - c.z*pm;
      r2[0] = r.x - r.z*pm;

      pm = a.y/a.z;
      m22[1][0] = b.y - b.z*pm;
      m22[1][1] = c.y - c.z*pm;
      r2[1] = r.y - r.z*pm;	      
      _Solve2x2(r_y, r_z, m22, r2);
      r_x = (r.z - b.z*r_y - c.z*r_z)/a.z;
    }
  } 
}
#endif	// tledCUDA_operators_CU
