#ifndef tledCUDA_operators_H
#define tledCUDA_operators_H

#ifdef _CUDA_5PLUS_SDK
#include <helper_math.h>
#endif

__device__ float3 operator+(const float3 a, const float3 b);
__device__ float3 operator+(const float3 a, const float b);
__device__ float3 operator+(const float3 &a, const float4 &b);
__device__ float3 operator+(const float a, const float3 b);
__device__ float2 operator+(const float2 &a, const float2 &b);
__device__ float3 operator+(const float3 a, const int b);
__device__ float3 operator+(const int a, const float3 b);
__device__ int3 operator+(const int3 a, const int b);
__device__ int3 operator+(const int a, const int3 b);
__device__ float4 operator+(const float4 a, const float4 b);
__device__ int4 operator+(const int4 a, const int b);

#ifndef _CUDA_5PLUS_SDK
__device__ void operator+=(float3 &r_a, const float3 &b);
__device__ void operator+=(float4 &a, const float4 &b);
__device__ void operator+=(float4 &a, const float3 &b);
__device__ void operator+=(float2 &a, const float2 &b);
#endif

__device__ float3 operator-(const float3 a, const float3 b);
__device__ float3 operator-(const float3 a, const float b);
__device__ float3 operator-(const float a, const float3 b);
__device__ float3 operator-(const float3 a, const int b);
__device__ float3 operator-(const int a, const float3 b);
__device__ float3 operator-(const float3 a);
__device__ float4 operator-(const float4 a, const float4 b);
__device__ float4 operator-(const float4 a);

#ifndef _CUDA_5PLUS_SDK
__device__ void operator-=(float3 &r_a, const float3 &b);
__device__ void operator-=(float4 &r_a, const float3 &b);
#endif

__device__ float3 operator*(const float a, const float3 b);
__device__ float3 operator*(const float3 a, const float b);
__device__ float3 operator*(const float3 a, const float3 b);
__device__ int3 operator*(const int3 a, const int3 b);
__device__ float3 operator*(const float3 a, const int3 b);
__device__ float3 operator*(const int3 a, const float3 b);
__device__ float4 operator*(const float a, const float4 b);
__device__ float4 operator*(const float4 a, const int4 b);
__device__ float4 operator*(const int4 a, const float4 b);

__device__ void operator*=(float3 &r_a, const float s);
__device__ void operator*=(float2 &r_a, const float s);

__device__ float4 operator/(const float4 a, const float b);
__device__ float3 operator/(const float3 a, const float b);

__device__ void operator/=(float3 &r_a, const float s);
__device__ void operator/=(float2 &r_a, const float s);
__device__ int3 operator%(const int3 a, const int b);

__device__ int4 operator%(const int4 a, const int b);

__device__ bool operator==(const int2 &a, const int2 &b);
__device__ bool operator!=(const int2 &a, const int2 &b);

__device__ float dot(const float3 &a, const float4 &b);
__device__ float dot(const float4 &a, const float3 &b);

#ifndef _CUDA_5PLUS_SDK
__device__ float dot(const float4 a, const float4 b);
__device__ float dot(const float3 a, const float3 b);
__device__ float dot(const float2 &a, const float2 &b);

__device__ float3 cross(const float3 &a, const float3 &b);
#endif

/** Interpolates between a, b based on a parameter t in [0, 1] (no checks performed on parameter!). */
__device__ float3 interpolate(const float3 &a, const float3 &b, const float t);

__device__ float norm(const float3 &a);
__device__ float norm(const float2 &a);

namespace tledCUDAMaths {
/**
 * \brief PLU direct solver.
 *
 *
 * Solves:<br>
 * \f$(a|b|c)\cdot(x,y,z)^T = r\f$
 */
  __device__ void SolveEquation3x3(float &r_x, float &r_y, float &r_z, const float3 &a, const float3 &b, const float3 &c, const float3 &r);
}
#endif
