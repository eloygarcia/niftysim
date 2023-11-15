// =========================================================================
// File:       tledSolverGPU_ROM_kernels.cu
// Purpose:    CUDA kernels for tledSolverGPU_ROM class
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   cu
// Created:    April 2011
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledSolverGPU_ROM_kernels_CU
#define tledSolverGPU_ROM_kernels_CU

#include "tledCUDA_operators.cu"
#include "tledContactManager.h"
#include "tledDeviceDeclarations.h"
#include "float.h"

namespace tledSolverGPU_ROM_kernels
{

// Compute 2nd Piola-Kirchhoff stress - linear elastic model
__device__ void SPK_LE(const unsigned int tid, float3* C_a, float3* C_b, float3* SPK_a, float3* SPK_b)
{
   float4 tempF4 = tex1Dfetch(txElasticParams,tid);	// M = tempF4.x, K = tempF4.y
   // 2nd Piola-Kirchhoff stress
   *SPK_a = make_float3(tempF4.x*(*C_a).x + tempF4.y*((*C_a).y + (*C_a).z) - tempF4.z,
                        tempF4.x*(*C_a).y + tempF4.y*((*C_a).x + (*C_a).z) - tempF4.z,
                        tempF4.x*(*C_a).z + tempF4.y*((*C_a).x + (*C_a).y) - tempF4.z);
   *SPK_b = make_float3(tempF4.w*(*C_b).x,
                        tempF4.w*(*C_b).y,
                        tempF4.w*(*C_b).z);
}

// Compute 2nd Piola-Kirchhoff stress - neo-Hookean hyperelastic model
__device__ void SPK_NH(const unsigned int tid, float3* C_a, float3* C_b, float J, float3* SPK_a, float3* SPK_b)
{
   float4 tempF4 = tex1Dfetch(txElasticParams,tid);	// M = tempF4.x, K = tempF4.y
   float x1 = powf(J,-2.0f/3.0f)*tempF4.x;
   float x2 = (tempF4.y*J*(J-1.0f) - x1*((*C_a).x+(*C_a).y+(*C_a).z)/3.0f)/
            ((*C_a).x*((*C_a).y*(*C_a).z - (*C_b).y*(*C_b).y) - (*C_a).y*(*C_b).z*(*C_b).z - 
            (*C_b).x*((*C_a).z*(*C_b).x - 2.0f*(*C_b).z*(*C_b).y));
   //Note: detC = C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23);

   // 2nd Piola-Kirchhoff stress
   *SPK_a = make_float3(((*C_a).y*(*C_a).z - (*C_b).y*(*C_b).y)*x2 + x1,
                           ((*C_a).x*(*C_a).z - (*C_b).z*(*C_b).z)*x2 + x1,
                           ((*C_a).x*(*C_a).y - (*C_b).x*(*C_b).x)*x2 + x1);
   *SPK_b = make_float3(((*C_b).z*(*C_b).y - (*C_b).x*(*C_a).z)*x2,
                           ((*C_b).z*(*C_b).x - (*C_b).y*(*C_a).x)*x2,
                           ((*C_b).x*(*C_b).y - (*C_b).z*(*C_a).y)*x2);
}

// Compute 2nd Piola-Kirchhoff stress - transversely isotropic hyperelastic model
__device__ void SPK_TI(const unsigned int tid, float3* C_a, float3* C_b, float J, float3* SPK_a, float3* SPK_b)
{
   // 2 MParams arrays -> need all of 1st, but only first two of 2nd. Get 2nd first, and store in a float2
   // variable, then get 1st -> confusing, but minimises registers
   float4 tempF4 = tex1Dfetch(txElasticParams,2*tid+1);
   float2 tempF2 = make_float2(tempF4.x,tempF4.y);
   tempF4 = tex1Dfetch(txElasticParams,2*tid);
   // Resulting param list:
   // Mu = tempF4.x
   // Kappa = tempF4.y
   // Eta = tempF4.z
   // a0 = tempF4.w
   // a1 = tempF2.x
   // a2 = tempF2.y

   float x5 = tempF4.y*J*(J-1);
   J = powf(J,-2.0f/3.0f);	// J = J^(-2/3)
   float x1 = J*tempF4.x;
   float x2 = J*(tempF4.w*tempF4.w*(*C_a).x+tempF2.x*tempF2.x*(*C_a).y+tempF2.y*tempF2.y*(*C_a).z
                  +2.0f*tempF4.w*tempF2.x*(*C_b).x+2.0f*tempF2.x*tempF2.y*(*C_b).y+2.0f*tempF4.w*tempF2.y*(*C_b).z) - 1.0f; // Bracketed term is I4 = A:C
   float x3 = J*tempF4.z*x2;
   float x4 = -(tempF4.z*x2*(x2+1) + x1*((*C_a).x+(*C_a).y+(*C_a).z))/3.0f;

   // C inverse
   float invdetC = 1.0f/((*C_a).x*((*C_a).y*(*C_a).z - (*C_b).y*(*C_b).y) - 
            (*C_a).y*(*C_b).z*(*C_b).z - (*C_b).x*((*C_a).z*(*C_b).x - 2.0f*(*C_b).z*(*C_b).y));

   // Total stresses
   // Si terms are isochoric, Sv terms are volumetric
   //SPK_a = make_float3(Si11 + Sv11,Si22 + Sv22,Si33 + Sv33);
   *SPK_a = make_float3(x3*tempF4.w*tempF4.w + (x4+x5)*((*C_a).y*(*C_a).z - (*C_b).y*(*C_b).y)*invdetC + x1,
                        x3*tempF2.x*tempF2.x + (x4+x5)*((*C_a).x*(*C_a).z - (*C_b).z*(*C_b).z)*invdetC + x1,
                        x3*tempF2.y*tempF2.y + (x4+x5)*((*C_a).x*(*C_a).y - (*C_b).x*(*C_b).x)*invdetC + x1);
   //SPK_b = make_float3(Si12 + Sv12,Si23 + Sv23,Si13 + Sv13);
   *SPK_b = make_float3(x3*tempF4.w*tempF2.x + (x4+x5)*((*C_b).z*(*C_b).y - (*C_b).x*(*C_a).z)*invdetC,
                        x3*tempF2.x*tempF2.y + (x4+x5)*((*C_b).x*(*C_b).z - (*C_b).y*(*C_a).x)*invdetC,
                        x3*tempF4.w*tempF2.y + (x4+x5)*((*C_b).x*(*C_b).y - (*C_b).z*(*C_a).y)*invdetC);
}

// Compute 2nd Piola-Kirchhoff stress - Arruda-Boyce hyperelastic model
__device__ void SPK_AB(const unsigned int tid, float3* C_a, float3* C_b, float J, float3* SPK_a, float3* SPK_b)
{
   float4 tempF4 = tex1Dfetch(txElasticParams,tid);	// M = tempF4.x, Lm = tempF4.y, K = tempF4.z
   float I1b = powf(J,-2.0f/3.0f)*((*C_a).x+(*C_a).y+(*C_a).z);
   // c1 = 1/2; c2 = 1/20; c3 = 11/1050; c4 = 19/7000; c5 = 519/673750;
   float Lm2 = tempF4.y*tempF4.y;
   float c1 = 0.5f;
   float c2 = 0.05f; c2 *= 2.0f*I1b/Lm2;
   float I1b2 = I1b*I1b;
   Lm2 *= tempF4.y*tempF4.y; // Lm2 = Lm^4
   float c3 = 0.01047619047619f; c3 *= 3.0f*I1b2/Lm2;
   I1b2 *= I1b; // I1b2 = I1b^3
   Lm2 *= tempF4.y*tempF4.y; // Lm2 = Lm^6
   float c4 = 0.002714285714286f; c4 *= 4.0f*I1b2/Lm2;
   I1b2 *= I1b; // I1b2 = I1b^4
   Lm2 *= tempF4.y*tempF4.y; // Lm2 = Lm^8
   float c5 = 7.703153988868275e-4f; c5 *= 5.0f*I1b2/Lm2;
   float x1 = powf(J,-2.0f/3.0f)*2*tempF4.x*(c1 + c2 + c3 + c4 + c5);
   float x2 = (0.5f*tempF4.z*J*(J - 1.0f/J) - x1*((*C_a).x+(*C_a).y+(*C_a).z)/3.0f)/
            ((*C_a).x*((*C_a).y*(*C_a).z - (*C_b).y*(*C_b).y) - (*C_a).y*(*C_b).z*(*C_b).z - 
            (*C_b).x*((*C_a).z*(*C_b).x - 2.0f*(*C_b).z*(*C_b).y));
   //Note: detC = C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23);

   // 2nd Piola-Kirchhoff stress
   *SPK_a = make_float3(((*C_a).y*(*C_a).z - (*C_b).y*(*C_b).y)*x2 + x1,
                           ((*C_a).x*(*C_a).z - (*C_b).z*(*C_b).z)*x2 + x1,
                           ((*C_a).x*(*C_a).y - (*C_b).x*(*C_b).x)*x2 + x1);
   *SPK_b = make_float3(((*C_b).z*(*C_b).y - (*C_b).x*(*C_a).z)*x2,
                           ((*C_b).z*(*C_b).x - (*C_b).y*(*C_a).x)*x2,
                           ((*C_b).x*(*C_b).y - (*C_b).z*(*C_a).y)*x2);

}

// Compute 2nd Piola-Kirchhoff stress - polynomial hyperelastic model
__device__ void SPK_PY(const unsigned int tid, float3* C_a, float3* C_b, float J, float3* SPK_a, float3* SPK_b)
{
   // 2 MParams arrays -> need all of 1st, but only first two of 2nd. Get 2nd first, and store in a float2
   // variable, then get 1st -> confusing, but minimises registers
   float4 tempF4 = tex1Dfetch(txElasticParams,2*tid+1);
   float2 tempF2 = make_float2(tempF4.x,tempF4.y);
   tempF4 = tex1Dfetch(txElasticParams,2*tid);
   // Resulting param list:
   // c10 = tempF4.x
   // c01 = tempF4.y
   // c20 = tempF4.z
   // c02 = tempF4.w
   // c11 = tempF2.x
   // Kappa = tempF2.y
   
   float J23 = powf(J,-2.0f/3.0f);
   // Determinant of C
   float detC = ((*C_a).x*((*C_a).y*(*C_a).z - (*C_b).y*(*C_b).y) - 
                  (*C_a).y*(*C_b).z*(*C_b).z - (*C_b).x*((*C_a).z*(*C_b).x
                  - 2.0f*(*C_b).z*(*C_b).y));
   // Invariants
   float I1 = (*C_a).x + (*C_a).y + (*C_a).z;
   float I2 = 0.5f*(I1*I1 - ((*C_a).x*(*C_a).x + 2.0f*(*C_b).x*(*C_b).x + 2.0f*(*C_b).z*(*C_b).z + 
                           (*C_a).y*(*C_a).y + 2.0f*(*C_b).y*(*C_b).y + (*C_a).z*(*C_a).z));
   float I1b = J23*I1;
   float I2b = J23*J23*I2;
   
   // Some convenience variables
   float a1 = I1b-3.0f;
   float a2 = I2b-3.0f;
   float a3 = 2.0f*J23*I2/3.0f;
   float T1 = 2.0f*J23*(tempF4.x + tempF4.y*I1b + 2.0f*tempF4.z*a1 + 2.0f*tempF4.w*a2*I1b 
                     + tempF2.x*(a2 + a1*I1b));
   float T2 = -2.0f*J23*J23*(tempF4.y + 2.0f*tempF4.w*a2 + tempF2.x*a1);
   float T3 = ( tempF2.y*J*(J-1.0f)
               -2.0f*J23*(tempF4.x*I1/3.0f + tempF4.y*a3 + 2.0f*tempF4.z*a1*I1/3.0f + 2.0f*tempF4.w*a2*a3 
               + tempF2.x*(a2*I1/3.0f + a1*a3))
               )/detC;

   // 2nd Piola-Kirchhoff stress
   *SPK_a = make_float3( (*C_a).x*T2 + ((*C_a).y*(*C_a).z - (*C_b).y*(*C_b).y)*T3 + T1,
                         (*C_a).y*T2 + ((*C_a).x*(*C_a).z - (*C_b).z*(*C_b).z)*T3 + T1,
                         (*C_a).z*T2 + ((*C_a).x*(*C_a).y - (*C_b).x*(*C_b).x)*T3 + T1 );
   *SPK_b = make_float3( (*C_b).x*T2 + ((*C_b).z*(*C_b).y - (*C_b).x*(*C_a).z)*T3,
                         (*C_b).y*T2 + ((*C_b).z*(*C_b).x - (*C_b).y*(*C_a).x)*T3,
                         (*C_b).z*T2 + ((*C_b).x*(*C_b).y - (*C_b).z*(*C_a).y)*T3 );
}

// Compute 2nd Piola-Kirchhoff stress - neo-Hookean visco-hyperelastic model
__device__ void SPK_NHV(const unsigned int tid, float3* C_a, float3* C_b, float J,
                        float4* g_StressStateIso, float4* g_StressStateVol, float3* SPK_a, float3* SPK_b)
{
   float C11 = (*C_a).x;
   float C22 = (*C_a).y;
   float C33 = (*C_a).z;
   float C12 = (*C_b).x;
   float C23 = (*C_b).y;
   float C13 = (*C_b).z;

   float4 tempF4 = tex1Dfetch(txElasticParams,tid);	// M = tempF4.x, K = tempF4.y
   float x1 = powf(J,-2.0f/3.0f)*tempF4.x;
   float invdetC = 1.0f/(C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2.0f*C13*C23));
   float x2 = -x1*(C11+C22+C33)/3.0f*invdetC;
   float x3 = tempF4.y*J*(J-1.0f)*invdetC;

   // Elastic components of response *********************************
   // SiE
   float SiE11 = x2*(C22*C33 - C23*C23) + x1;
   float SiE22 = x2*(C11*C33 - C13*C13) + x1;
   float SiE33 = x2*(C11*C22 - C12*C12) + x1;
   float SiE12 = x2*(C13*C23 - C12*C33);
   float SiE23 = x2*(C12*C13 - C11*C23);
   float SiE13 = x2*(C12*C23 - C13*C22);
   // SvE
   float SvE11 = x3*(C22*C33 - C23*C23);
   float SvE22 = x3*(C11*C33 - C13*C13);
   float SvE33 = x3*(C11*C22 - C12*C12);
   float SvE12 = x3*(C13*C23 - C12*C33);
   float SvE23 = x3*(C12*C13 - C11*C23);
   float SvE13 = x3*(C12*C23 - C13*C22);

   // State variables ************************************************
   // Reuse C-variables for summing state variables
   // --> not required anymore, and saves declaring another 6 floats for this purpose
   C11 = 0.0f; C22 = 0.0f; C33 = 0.0f; C12 = 0.0f; C23 = 0.0f; C13 = 0.0f;
   int2 NumProny = tex1Dfetch(txNumProny,tid);
   float2 euler;
   // Isochoric
   for (int i = 0; i < NumProny.x; i++)
   {
      tempF4 = tex1Dfetch(txStateIso,tid*2*c_maxNumViscTerms.x + i);	// Get state variables
      euler = tex1Dfetch(txEulerIso,tid*c_maxNumViscTerms.x + i);	// Get backward Euler coeffs
      tempF4.x *= euler.y; tempF4.x += euler.x*SiE11; C11 += tempF4.x;
      tempF4.y *= euler.y; tempF4.y += euler.x*SiE22; C22 += tempF4.y;
      tempF4.z *= euler.y; tempF4.z += euler.x*SiE33; C33 += tempF4.z;
      tempF4.w *= euler.y; tempF4.w += euler.x*SiE12; C12 += tempF4.w;
      g_StressStateIso[tid*2*c_maxNumViscTerms.x + i] = tempF4;	// Write updated state variables to memory
      
      tempF4 = tex1Dfetch(txStateIso,tid*2*c_maxNumViscTerms.x + i + 1);	// Get state variables
      tempF4.x *= euler.y; tempF4.x += euler.x*SiE23; C23 += tempF4.x;
      tempF4.y *= euler.y; tempF4.y += euler.x*SiE13; C13 += tempF4.y;
      tempF4.z = 0.0f;
      tempF4.w = 0.0f;
      g_StressStateIso[tid*2*c_maxNumViscTerms.x + i + 1] = tempF4;
   }
   // Volumetric
   for (int i = 0; i < NumProny.y; i++)
   {
      tempF4 = tex1Dfetch(txStateVol,tid*2*c_maxNumViscTerms.y + i);
      euler = tex1Dfetch(txEulerVol,tid*c_maxNumViscTerms.y + i);
      tempF4.x *= euler.y; tempF4.x += euler.x*SvE11; C11 += tempF4.x;
      tempF4.y *= euler.y; tempF4.y += euler.x*SvE22; C22 += tempF4.y;
      tempF4.z *= euler.y; tempF4.z += euler.x*SvE33; C33 += tempF4.z;
      tempF4.w *= euler.y; tempF4.w += euler.x*SvE12; C12 += tempF4.w;
      g_StressStateVol[tid*2*c_maxNumViscTerms.y + i] = tempF4;
      
      tempF4 = tex1Dfetch(txStateVol,tid*2*c_maxNumViscTerms.y + i + 1);
      tempF4.x *= euler.y; tempF4.x += euler.x*SvE23; C23 += tempF4.x;
      tempF4.y *= euler.y; tempF4.y += euler.x*SvE13; C13 += tempF4.y;
      tempF4.z = 0.0f;
      tempF4.w = 0.0f;
      g_StressStateVol[tid*2*c_maxNumViscTerms.y + i + 1] = tempF4;
   }

   // Total stress ***************************************************
   *SPK_a = make_float3(SiE11 + SvE11 - C11,
                        SiE22 + SvE22 - C22,
                        SiE33 + SvE33 - C33);
   *SPK_b = make_float3(SiE12 + SvE12 - C12,
                        SiE23 + SvE23 - C23,
                        SiE13 + SvE13 - C13);
}

// Compute 2nd Piola-Kirchhoff stress - transversely isotropic visco-hyperelastic model
__device__ void SPK_TIV(const unsigned int tid, float3* C_a, float3* C_b, float J,
                        float4* g_StressStateIso, float4* g_StressStateVol, float3* SPK_a, float3* SPK_b)
{
   float C11 = (*C_a).x;
   float C22 = (*C_a).y;
   float C33 = (*C_a).z;
   float C12 = (*C_b).x;
   float C23 = (*C_b).y;
   float C13 = (*C_b).z;

   // 2 MParams arrays -> need all of 1st, but only first two of 2nd. Get 2nd first, and store in a float2
   // variable, then get 1st -> confusing, but minimises registers
   float4 tempF4 = tex1Dfetch(txElasticParams,2*tid+1);
   float2 tempF2 = make_float2(tempF4.x,tempF4.y);
   tempF4 = tex1Dfetch(txElasticParams,2*tid);
   // Resulting param list:
   // Mu = tempF4.x
   // Kappa = tempF4.y
   // Eta = tempF4.z
   // a0 = tempF4.w
   // a1 = tempF2.x
   // a2 = tempF2.y

   float invdetC = 1.0f/(C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2.0f*C13*C23));
   float x5 = ( tempF4.y*J*(J-1.0f) )*invdetC;
   J = powf(J,-2.0f/3.0f);	// J = J^(-2/3)
   float x1 = J*tempF4.x;
   float x2 = J*(tempF4.w*tempF4.w*C11+tempF2.x*tempF2.x*C22+tempF2.y*tempF2.y*C33
            +2.0f*tempF4.w*tempF2.x*C12+2.0f*tempF2.x*tempF2.y*C23+2.0f*tempF4.w*tempF2.y*C13) - 1; // Bracketed term is I4 = A:C
   float x3 = J*tempF4.z*x2;
   float x4 = ( -(tempF4.z*x2*(x2+1.0f) + x1*(C11+C22+C33))/3.0f )*invdetC;

   // Elastic components of response *********************************
   // SiE
   float SiE11 = x3*tempF4.w*tempF4.w + x4*(C22*C33 - C23*C23) + x1;
   float SiE22 = x3*tempF2.x*tempF2.x + x4*(C11*C33 - C13*C13) + x1;
   float SiE33 = x3*tempF2.y*tempF2.y + x4*(C11*C22 - C12*C12) + x1;
   float SiE12 = x3*tempF4.w*tempF2.x + x4*(C13*C23 - C12*C33);
   float SiE23 = x3*tempF2.x*tempF2.y + x4*(C12*C13 - C11*C23);
   float SiE13 = x3*tempF4.w*tempF2.y + x4*(C12*C23 - C13*C22);
   // SvE
   float SvE11 = x5*(C22*C33 - C23*C23);
   float SvE22 = x5*(C11*C33 - C13*C13);
   float SvE33 = x5*(C11*C22 - C12*C12);
   float SvE12 = x5*(C13*C23 - C12*C33);
   float SvE23 = x5*(C12*C13 - C11*C23);
   float SvE13 = x5*(C12*C23 - C13*C22);

   // State variables ************************************************
   // Reuse C-variables for summing state variables
   // --> not required anymore, and saves declaring another 6 floats for this purpose
   C11 = 0.0f; C22 = 0.0f; C33 = 0.0f; C12 = 0.0f; C23 = 0.0f; C13 = 0.0f;
   int2 NumProny = tex1Dfetch(txNumProny,tid);
   // Isochoric
   for (int i = 0; i < NumProny.x; i++)
   {
      tempF4 = tex1Dfetch(txStateIso,tid*2*c_maxNumViscTerms.x + i);	// Get state variables
      tempF2 = tex1Dfetch(txEulerIso,tid*c_maxNumViscTerms.x + i);	// Get backward Euler coeffs
      tempF4.x *= tempF2.y; tempF4.x += tempF2.x*SiE11; C11 += tempF4.x;
      tempF4.y *= tempF2.y; tempF4.y += tempF2.x*SiE22; C22 += tempF4.y;
      tempF4.z *= tempF2.y; tempF4.z += tempF2.x*SiE33; C33 += tempF4.z;
      tempF4.w *= tempF2.y; tempF4.w += tempF2.x*SiE12; C12 += tempF4.w;
      g_StressStateIso[tid*2*c_maxNumViscTerms.x + i] = tempF4;	// Write updated state variables to memory
      
      tempF4 = tex1Dfetch(txStateIso,tid*2*c_maxNumViscTerms.x + i + 1);	// Get state variables
      tempF4.x *= tempF2.y; tempF4.x += tempF2.x*SiE23; C23 += tempF4.x;
      tempF4.y *= tempF2.y; tempF4.y += tempF2.x*SiE13; C13 += tempF4.y;
      tempF4.z = 0.0f;
      tempF4.w = 0.0f;
      g_StressStateIso[tid*2*c_maxNumViscTerms.x + i + 1] = tempF4;
   }
   // Volumetric
   for (int i = 0; i < NumProny.y; i++)
   {
      tempF4 = tex1Dfetch(txStateVol,tid*2*c_maxNumViscTerms.y + i);
      tempF2 = tex1Dfetch(txEulerVol,tid*c_maxNumViscTerms.y + i);
      tempF4.x *= tempF2.y; tempF4.x += tempF2.x*SvE11; C11 += tempF4.x;
      tempF4.y *= tempF2.y; tempF4.y += tempF2.x*SvE22; C22 += tempF4.y;
      tempF4.z *= tempF2.y; tempF4.z += tempF2.x*SvE33; C33 += tempF4.z;
      tempF4.w *= tempF2.y; tempF4.w += tempF2.x*SvE12; C12 += tempF4.w;
      g_StressStateVol[tid*2*c_maxNumViscTerms.y + i] = tempF4;
      
      tempF4 = tex1Dfetch(txStateVol,tid*2*c_maxNumViscTerms.y + i + 1);
      tempF4.x *= tempF2.y; tempF4.x += tempF2.x*SvE23; C23 += tempF4.x;
      tempF4.y *= tempF2.y; tempF4.y += tempF2.x*SvE13; C13 += tempF4.y;
      tempF4.z = 0.0f;
      tempF4.w = 0.0f;
      g_StressStateVol[tid*2*c_maxNumViscTerms.y + i + 1] = tempF4;
   }

   // Total stress ***************************************************
   *SPK_a = make_float3(SiE11 + SvE11 - C11,
                        SiE22 + SvE22 - C22,
                        SiE33 + SvE33 - C33);
   *SPK_b = make_float3(SiE12 + SvE12 - C12,
                        SiE23 + SvE23 - C23,
                        SiE13 + SvE13 - C13);
}

// Compute the barycentric coordinates of the point P in the triangle T
__device__ void TriangleBaryCoords(float3* T1, float3* T2, float3* T3, float3* p, float3* uvw)
{
   // v0 = T2 - T1;
   // v1 = T3 - T1;
   // v2 = P - T1;
   float d00 = dot((*T2)-(*T1),(*T2)-(*T1));
   float d01 = dot((*T2)-(*T1),(*T3)-(*T1));
   float d11 = dot((*T3)-(*T1),(*T3)-(*T1));
   float d20 = dot((*p)-(*T1),(*T2)-(*T1));
   float d21 = dot((*p)-(*T1),(*T3)-(*T1));
   float denom = 1/(d00*d11 - d01*d01);
   uvw->y = (d11*d20 - d01*d21)*denom;
   uvw->z = (d00*d21 - d01*d20)*denom;
   uvw->x = 1 - uvw->y - uvw->z;
}

// Check for penetration of contact cylinder master surface by slave nodes
__device__ void ResolveCylinderContact(const float4* p, const float4* cyl_orig, const float4* cyl_axis, float3* Unext_in)
{
   // NB: cylinder params arranged as
   // orig = (cyl_orig.x,cyl_orig.y,cyl_orig.z)
   // radius = cyl_orig.w
   // axis = (cyl_axis.x,cyl_axis.y,cyl_axis.z)
   // length = cyl_axis.w
   
   // Transform p into cylinder coords
   float3 P = make_float3((*p).x,(*p).y,(*p).z) - make_float3((*cyl_orig).x,(*cyl_orig).y,(*cyl_orig).z); // P = p - orig
   // Project onto axis
   float nb = dot(P,make_float3((*cyl_axis).x,(*cyl_axis).y,(*cyl_axis).z)); // nb = dot(P,axis)
   if ( (nb <= 0) | (nb >= (*cyl_axis).w) ) // Just return if p is without axis limits
      return;
   float3 b = nb*make_float3((*cyl_axis).x,(*cyl_axis).y,(*cyl_axis).z); // b = nb*axis
   // Compute radial distance to p
   float3 B = P - b; // Vector from axis to P
   float nB = dot(B,B); // Squared distance from axis to P
   if (nB < (*cyl_orig).w * (*cyl_orig).w) // Check if within cylinder surface
   {
      // Closest point: c = b + cyl_R*B/sqrt(nB)
      // Transform c to global coords: c = c + *cyl_orig
      // Vector from current to required pos: u = c - *p
      // Update Unext: Unext = Unext + u
      // Combine into 1 statement:
      nB = sqrt(nB);
      *Unext_in = *Unext_in + b + (*cyl_orig).w*B/nB + make_float3((*cyl_orig).x,(*cyl_orig).y,(*cyl_orig).z) - make_float3((*p).x,(*p).y,(*p).z);
   }
}

// Check for penetration of contact ultrasound probe master surface by slave nodes
__device__ void ResolveUSProbeContact(const float4* p, const float4* prb_orig, const float4* prb_axis, float3* Unext_in)
{
   // NB: probe params arranged as
   // orig = (prb_orig.x,prb_orig.y,prb_orig.z)
   // radius = prb_orig.w
   // axis = (prb_axis.x,prb_axis.y,prb_axis.z)
   // length = prb_axis.w
      
   // Transform p into cylinder coords
   float3 P = make_float3((*p).x,(*p).y,(*p).z) - make_float3((*prb_orig).x,(*prb_orig).y,(*prb_orig).z); // P = p - orig
   // Project onto axis
   float nb = dot(P,make_float3((*prb_axis).x,(*prb_axis).y,(*prb_axis).z)); // nb = dot(P,axis)
   if (nb <= 0) // Just return if p is below probe bottom
      return;
   else if ( (nb >= (*prb_axis).w) ) // If p is beyond cylinder top, check if it is within the spherical tip
   {
      // Transform p into sphere coords
      P = P - (*prb_axis).w * make_float3((*prb_axis).x,(*prb_axis).y,(*prb_axis).z); // Back to global: P=P+orig; Into sphere: P=P-(orig+L*axis); Thus: P=P+orig-orig-L*axis=P-L*axis
      float nP = dot(P,P); // Squared distance from sphere centre to P
      if (nP >= (*prb_orig).w * (*prb_orig).w) // Check if P is within sphere surface
         return;
      
      // Closest point: c = R*P/sqrt(nP)
      // Transform c into global coords: c = c + orig L*axis = R*P/sqrt(nP) + orig + L*axis
      // Vector from current to req'd pos: u = c - p = R*P/sqrt(nP) + orig + L*axis - p
      // Update Unext: Unext_in = Unext_in + u
      //                        = Unext_in + R*P/sqrt(nP) + orig + L*axis - p
      nP = sqrt(nP);
      *Unext_in = *Unext_in + (*prb_orig).w*P/nP + make_float3((*prb_orig).x,(*prb_orig).y,(*prb_orig).z) 
                  + (*prb_axis).w * make_float3((*prb_axis).x,(*prb_axis).y,(*prb_axis).z)
                  - make_float3((*p).x,(*p).y,(*p).z);
   }
   else // p is within limits of cylinder
   {
      float3 b = nb * make_float3((*prb_axis).x,(*prb_axis).y,(*prb_axis).z);
      // Compute radial distance to p
      float3 B = P-b; // Vector from axis to P
      float nB = dot(B,B); // Squared
      if (nB < (*prb_orig).w * (*prb_orig).w) // Check if within cylinder surface
      {
         // Closest point: c = b + cyl_R*B/sqrt(nB)
         // Transform c to global coords: c = c + *cyl_orig
         // Vector from current to required pos: u = c - *p
         // Update Unext: Unext = Unext + u
         // Combine into 1 statement:
         nB = sqrt(nB);
         *Unext_in = *Unext_in + b + (*prb_orig).w*B/nB + make_float3((*prb_orig).x,(*prb_orig).y,(*prb_orig).z) - make_float3((*p).x,(*p).y,(*p).z);
      }
   }
}

// Check for penetration of contact plate master surface by slave nodes
__device__ void ResolvePlateContact(const float4* p, const float4* a, const float4* b, const float4* c, float3* Unext_in)
{
   // Surface norm
   float3 ab = make_float3((*b).x-(*a).x,(*b).y-(*a).y,(*b).z-(*a).z);
   float3 ac = make_float3((*c).x-(*a).x,(*c).y-(*a).y,(*c).z-(*a).z);
   float3 N = make_float3(ab.y*ac.z - ab.z*ac.y,
                           ab.z*ac.x - ab.x*ac.z,
                           ab.x*ac.y - ab.y*ac.x); // N = cross(ab,ac)
   N = N/sqrt(dot(N,N)); // Normalise
   // Shift p to local coords
   float3 d = make_float3((*p).x,(*p).y,(*p).z)-make_float3((*a).x,(*a).y,(*a).z);
   // Signed distance from node to plane of the plate
   float D = dot(N,d);
   if (D >= 0) // No penetration
      return;
   // Project onto local axes
   float L1 = sqrt(dot(ab,ab));
   if (dot(d,ab/L1) > L1) // Outside x-limit
      return;
   float L2 = sqrt(dot(ac,ac));
   if (dot(d,ac/L2) > L2) // Outside y-limit
      return;
   // Node has penetrated plate --> find closest point C on the plate
   // In this case, this is just the projection of p onto the plane of the plate
   // C = p - D*N
   // The displacement from this point to the current node location is
   // u = C - p = -D*N
   // Then the required displacement for this step is
   // Unext = Unext + u = Unext - D*N
   *Unext_in = *Unext_in - D*N;
}

// Compute internal forces for the current configuration - linear tetrahedron
__global__ void ComputeNewForcesT4_kernel(float4* g_FEl, float4* g_StressStateIso, float4* g_StressStateVol,
                                          float4* g_SPKa, float4* g_SPKb, float4* g_Ca, float4* g_Cb)
{
	// Note on nomenclature:
	// Variables preceded by g_ reside in global memory

	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
if (tid < c_NumEls)
{
	float4 Vol_MType_K = tex1Dfetch(txVol_MType_K,tid);
	// Shape function derivatives
	float Dh[4][3];
	float4 tempF4;
	tempF4 = tex1Dfetch(txDhDx,3*tid);
	Dh[0][0] = tempF4.x; Dh[1][0] = tempF4.y; Dh[2][0] = tempF4.z; Dh[3][0] = tempF4.w;
	tempF4 = tex1Dfetch(txDhDx,3*tid + 1);
	Dh[0][1] = tempF4.x; Dh[1][1] = tempF4.y; Dh[2][1] = tempF4.z; Dh[3][1] = tempF4.w;
	tempF4 = tex1Dfetch(txDhDx,3*tid + 2);
	Dh[0][2] = tempF4.x; Dh[1][2] = tempF4.y; Dh[2][2] = tempF4.z; Dh[3][2] = tempF4.w;

	// Displacements
	float u[4][3];
	int4 EInd;
	EInd = tex1Dfetch(txEInd,tid);
   tempF4 = tex1Dfetch(txUcurr,EInd.x);
	u[0][0] = tempF4.x; u[0][1] = tempF4.y; u[0][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.y);
	u[1][0] = tempF4.x; u[1][1] = tempF4.y; u[1][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.z);
	u[2][0] = tempF4.x; u[2][1] = tempF4.y; u[2][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.w);
	u[3][0] = tempF4.x; u[3][1] = tempF4.y; u[3][2] = tempF4.z;

	// Deformation gradient (transpose of), X = Du + I, Du = u^T*Dh
	float XT00,XT01,XT02,XT10,XT11,XT12,XT20,XT21,XT22;
	XT00 = u[0][0]*Dh[0][0]+u[1][0]*Dh[1][0]+u[2][0]*Dh[2][0]+u[3][0]*Dh[3][0] + 1;
	XT01 = u[0][1]*Dh[0][0]+u[1][1]*Dh[1][0]+u[2][1]*Dh[2][0]+u[3][1]*Dh[3][0];
	XT02 = u[0][2]*Dh[0][0]+u[1][2]*Dh[1][0]+u[2][2]*Dh[2][0]+u[3][2]*Dh[3][0];
	XT10 = u[0][0]*Dh[0][1]+u[1][0]*Dh[1][1]+u[2][0]*Dh[2][1]+u[3][0]*Dh[3][1];
	XT11 = u[0][1]*Dh[0][1]+u[1][1]*Dh[1][1]+u[2][1]*Dh[2][1]+u[3][1]*Dh[3][1] + 1;
	XT12 = u[0][2]*Dh[0][1]+u[1][2]*Dh[1][1]+u[2][2]*Dh[2][1]+u[3][2]*Dh[3][1];
	XT20 = u[0][0]*Dh[0][2]+u[1][0]*Dh[1][2]+u[2][0]*Dh[2][2]+u[3][0]*Dh[3][2];
	XT21 = u[0][1]*Dh[0][2]+u[1][1]*Dh[1][2]+u[2][1]*Dh[2][2]+u[3][1]*Dh[3][2];
	XT22 = u[0][2]*Dh[0][2]+u[1][2]*Dh[1][2]+u[2][2]*Dh[2][2]+u[3][2]*Dh[3][2] + 1;

	// Stress
	float3 SPKa, SPKb;
	// Right Cauchy-Green deformation tensor
	float3 Ca = make_float3(XT00*XT00 + XT01*XT01 + XT02*XT02,   // C11
                                 XT10*XT10 + XT11*XT11 + XT12*XT12,   // C22
                                 XT20*XT20 + XT21*XT21 + XT22*XT22);  // C33
	float3 Cb = make_float3(XT00*XT10 + XT01*XT11 + XT02*XT12,   // C12
                                 XT10*XT20 + XT11*XT21 + XT12*XT22,   // C23
                                 XT00*XT20 + XT01*XT21 + XT02*XT22);  // C13
        
   // Save the deformation
   g_Ca[tid] = make_float4(Ca.x,Ca.y,Ca.z,0.0f);
   g_Cb[tid] = make_float4(Cb.x,Cb.y,Cb.z,0.0f);

	if (Vol_MType_K.y == 1)	// LE
	{
      tledSolverGPU_ROM_kernels::SPK_LE(tid,&Ca,&Cb,&SPKa,&SPKb);
	}
	else if (Vol_MType_K.y == 2)	// NH
	{
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NH(tid,&Ca,&Cb,J,&SPKa,&SPKb);
	}
	else if (Vol_MType_K.y == 3)	// TI
	{
		float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TI(tid,&Ca,&Cb,J,&SPKa,&SPKb);
	}
	else if (Vol_MType_K.y == 4)	// NHV
	{
		float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NHV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
	}
	else if (Vol_MType_K.y == 5)	// TIV
	{
		float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TIV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
	}
   else if (Vol_MType_K.y == 6)	// AB
	{
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_AB(tid,&Ca,&Cb,J,&SPKa,&SPKb);
	}
   else if (Vol_MType_K.y == 7)	// PY
	{
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_PY(tid,&Ca,&Cb,J,&SPKa,&SPKb);
	}
        
   // Save the stresses
   g_SPKa[tid] = make_float4(SPKa.x,SPKa.y,SPKa.z,0.0f);
   g_SPKb[tid] = make_float4(SPKb.x,SPKb.y,SPKb.z,0.0f);

	// Compute nodal forces: FEl = V*X*S*Dh'
   // Step 1: S = S*V
   SPKa = SPKa*Vol_MType_K.x;
   SPKb = SPKb*Vol_MType_K.x;
   // Step 2: A = X*S
   float A00 = XT00*SPKa.x+XT10*SPKb.x+XT20*SPKb.z;
   float A01 = XT00*SPKb.x+XT10*SPKa.y+XT20*SPKb.y;
   float A02 = XT00*SPKb.z+XT10*SPKb.y+XT20*SPKa.z;
   float A10 = XT01*SPKa.x+XT11*SPKb.x+XT21*SPKb.z;
   float A11 = XT01*SPKb.x+XT11*SPKa.y+XT21*SPKb.y;
   float A12 = XT01*SPKb.z+XT11*SPKb.y+XT21*SPKa.z;
   float A20 = XT02*SPKa.x+XT12*SPKb.x+XT22*SPKb.z;
   float A21 = XT02*SPKb.x+XT12*SPKa.y+XT22*SPKb.y;
   float A22 = XT02*SPKb.z+XT12*SPKb.y+XT22*SPKa.z;
   // Step 3: FEl = A*Dh'
   for (int i = 0; i < 4; i++)
   {
      float Dhi0 = Dh[i][0]; float Dhi1 = Dh[i][1]; float Dhi2 = Dh[i][2];
      tempF4.x = A00*Dhi0+A01*Dhi1+A02*Dhi2;
      tempF4.y = A10*Dhi0+A11*Dhi1+A12*Dhi2;
      tempF4.z = A20*Dhi0+A21*Dhi1+A22*Dhi2;

      g_FEl[tid*4 + i] = tempF4;
   }
} // tid < c_NumEls
}

// Compute element pressures for the current configuration - nodal-averaged pressure tetrahedron
__global__ void ComputeElementPressuresT4ANP_kernel(float4* g_FEl)
{
   // Note on nomenclature:
   // Variables preceded by g_ reside in global memory

   const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
if (tid < c_NumEls)
{
   float4 Vol_MType_K = tex1Dfetch(txVol_MType_K,tid);	// Vol = Vol_MType_K.x, K = Vol_MType_K.z
   // Shape function derivatives
   float Dh[4][3];
   float4 tempF4;
   tempF4 = tex1Dfetch(txDhDx,3*tid);
   Dh[0][0] = tempF4.x; Dh[1][0] = tempF4.y; Dh[2][0] = tempF4.z; Dh[3][0] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,3*tid + 1);
   Dh[0][1] = tempF4.x; Dh[1][1] = tempF4.y; Dh[2][1] = tempF4.z; Dh[3][1] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,3*tid + 2);
   Dh[0][2] = tempF4.x; Dh[1][2] = tempF4.y; Dh[2][2] = tempF4.z; Dh[3][2] = tempF4.w;

   // Displacements
   float u[4][3];
   int4 EInd;
   EInd = tex1Dfetch(txEInd,tid);
   tempF4 = tex1Dfetch(txUcurr,EInd.x);
   u[0][0] = tempF4.x; u[0][1] = tempF4.y; u[0][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.y);
   u[1][0] = tempF4.x; u[1][1] = tempF4.y; u[1][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.z);
   u[2][0] = tempF4.x; u[2][1] = tempF4.y; u[2][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.w);
   u[3][0] = tempF4.x; u[3][1] = tempF4.y; u[3][2] = tempF4.z;

   // Deformation gradient (transpose of), X = Du + I, Du = u^T*Dh
   float XT00,XT01,XT02,XT10,XT11,XT12,XT20,XT21,XT22;
   XT00 = u[0][0]*Dh[0][0]+u[1][0]*Dh[1][0]+u[2][0]*Dh[2][0]+u[3][0]*Dh[3][0] + 1;
   XT01 = u[0][1]*Dh[0][0]+u[1][1]*Dh[1][0]+u[2][1]*Dh[2][0]+u[3][1]*Dh[3][0];
   XT02 = u[0][2]*Dh[0][0]+u[1][2]*Dh[1][0]+u[2][2]*Dh[2][0]+u[3][2]*Dh[3][0];
   XT10 = u[0][0]*Dh[0][1]+u[1][0]*Dh[1][1]+u[2][0]*Dh[2][1]+u[3][0]*Dh[3][1];
   XT11 = u[0][1]*Dh[0][1]+u[1][1]*Dh[1][1]+u[2][1]*Dh[2][1]+u[3][1]*Dh[3][1] + 1;
   XT12 = u[0][2]*Dh[0][1]+u[1][2]*Dh[1][1]+u[2][2]*Dh[2][1]+u[3][2]*Dh[3][1];
   XT20 = u[0][0]*Dh[0][2]+u[1][0]*Dh[1][2]+u[2][0]*Dh[2][2]+u[3][0]*Dh[3][2];
   XT21 = u[0][1]*Dh[0][2]+u[1][1]*Dh[1][2]+u[2][1]*Dh[2][2]+u[3][1]*Dh[3][2];
   XT22 = u[0][2]*Dh[0][2]+u[1][2]*Dh[1][2]+u[2][2]*Dh[2][2]+u[3][2]*Dh[3][2] + 1;

   float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);

   // El pressure
   // Store P*Vol: P = kappa*(J-1), Vol = Vol_MType.x
   g_FEl[tid*4].w = Vol_MType_K.z*(J-1)*Vol_MType_K.x;
} // tid < c_NumEls
}

// Compute nodal-averaged pressures for the current configuration - nodal-averaged pressure tetrahedron
__global__ void ComputeNodalPressuresT4ANP_kernel(float* g_Pa)
{
   const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
if (tid < c_NumNodes)
{
   float Va = tex1Dfetch(txVa,tid);
   float4 tempF4;
   // Gather element pressures
   float Pa = 0.0f;
   int2 NodeMap = tex1Dfetch(txNodeMap,tid);
   int2 FCds;
   for (int i = 0; i < NodeMap.y; i++)
   {
      FCds = tex1Dfetch(txFCds,NodeMap.x + i);	// First val contains element number (all that is required)
      tempF4 = tex1Dfetch(txFEl,4*FCds.x);	// El pressures are stored in the .w val of first FEl for each El
      Pa += tempF4.w;
   }
   // Store Pa/Va
   g_Pa[tid] = Pa/Va;
} // tid < c_NumNodes
}

// Compute modified element forces for the current configuration - nodal-averaged pressure tetrahedron
__global__ void ComputeModifiedElementForcesT4ANP_kernel(float4* g_FEl, float4* g_StressStateIso, float4* g_StressStateVol,
                                                         float4* g_SPKa, float4* g_SPKb, float4* g_Ca, float4* g_Cb)
{
   const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
if (tid < c_NumEls)
{
   float4 Vol_MType_K = tex1Dfetch(txVol_MType_K,tid);	// Vol = Vol_MType_K.x, MType = Vol_MType_K.y, K = Vol_MType_K.z
   // Shape function derivatives
   float Dh[4][3];
   float4 tempF4;
   tempF4 = tex1Dfetch(txDhDx,3*tid);
   Dh[0][0] = tempF4.x; Dh[1][0] = tempF4.y; Dh[2][0] = tempF4.z; Dh[3][0] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,3*tid + 1);
   Dh[0][1] = tempF4.x; Dh[1][1] = tempF4.y; Dh[2][1] = tempF4.z; Dh[3][1] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,3*tid + 2);
   Dh[0][2] = tempF4.x; Dh[1][2] = tempF4.y; Dh[2][2] = tempF4.z; Dh[3][2] = tempF4.w;

   // Displacements
   float u[4][3];
   int4 EInd;
   EInd = tex1Dfetch(txEInd,tid);
   tempF4 = tex1Dfetch(txUcurr,EInd.x);
   u[0][0] = tempF4.x; u[0][1] = tempF4.y; u[0][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.y);
   u[1][0] = tempF4.x; u[1][1] = tempF4.y; u[1][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.z);
   u[2][0] = tempF4.x; u[2][1] = tempF4.y; u[2][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.w);
   u[3][0] = tempF4.x; u[3][1] = tempF4.y; u[3][2] = tempF4.z;

   // Deformation gradient (transpose of), X = Du + I, Du = u^T*Dh
   float XT00,XT01,XT02,XT10,XT11,XT12,XT20,XT21,XT22;
   XT00 = u[0][0]*Dh[0][0]+u[1][0]*Dh[1][0]+u[2][0]*Dh[2][0]+u[3][0]*Dh[3][0] + 1;
   XT01 = u[0][1]*Dh[0][0]+u[1][1]*Dh[1][0]+u[2][1]*Dh[2][0]+u[3][1]*Dh[3][0];
   XT02 = u[0][2]*Dh[0][0]+u[1][2]*Dh[1][0]+u[2][2]*Dh[2][0]+u[3][2]*Dh[3][0];
   XT10 = u[0][0]*Dh[0][1]+u[1][0]*Dh[1][1]+u[2][0]*Dh[2][1]+u[3][0]*Dh[3][1];
   XT11 = u[0][1]*Dh[0][1]+u[1][1]*Dh[1][1]+u[2][1]*Dh[2][1]+u[3][1]*Dh[3][1] + 1;
   XT12 = u[0][2]*Dh[0][1]+u[1][2]*Dh[1][1]+u[2][2]*Dh[2][1]+u[3][2]*Dh[3][1];
   XT20 = u[0][0]*Dh[0][2]+u[1][0]*Dh[1][2]+u[2][0]*Dh[2][2]+u[3][0]*Dh[3][2];
   XT21 = u[0][1]*Dh[0][2]+u[1][1]*Dh[1][2]+u[2][1]*Dh[2][2]+u[3][1]*Dh[3][2];
   XT22 = u[0][2]*Dh[0][2]+u[1][2]*Dh[1][2]+u[2][2]*Dh[2][2]+u[3][2]*Dh[3][2] + 1;

   float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
   J = powf(J,-1.0f/3.0f); // J = J^-1/3

   // Compute modified pressure
   float Pm = 0.0f;
   Pm += tex1Dfetch(txPa,EInd.x);
   Pm += tex1Dfetch(txPa,EInd.y);
   Pm += tex1Dfetch(txPa,EInd.z);
   Pm += tex1Dfetch(txPa,EInd.w);
   Pm /= 16;

   // Compute modified Jacobian
   float Jm = Pm/Vol_MType_K.z + 1; // Jm = Pm/K + 1
   // Compute modified def. grad.
   Jm = powf(Jm,1.0f/3.0f); // Jm = Jm^1/3
   Jm = Jm*J;	// Jm = Jm^1/3 * J^-1/3
   XT00 *= Jm;
   XT01 *= Jm;
   XT02 *= Jm;
   XT10 *= Jm;
   XT11 *= Jm;
   XT12 *= Jm;
   XT20 *= Jm;
   XT21 *= Jm;
   XT22 *= Jm;

   // Stress
   float3 SPKa, SPKb;
   // Right Cauchy-Green deformation tensor
   float3 Ca = make_float3(XT00*XT00 + XT01*XT01 + XT02*XT02,	// C11
                           XT10*XT10 + XT11*XT11 + XT12*XT12,	// C22
                           XT20*XT20 + XT21*XT21 + XT22*XT22);	// C33
   float3 Cb = make_float3(XT00*XT10 + XT01*XT11 + XT02*XT12,	// C12
                           XT10*XT20 + XT11*XT21 + XT12*XT22,	// C23
                           XT00*XT20 + XT01*XT21 + XT02*XT22);	// C13
        
   // Save the deformation
   g_Ca[tid] = make_float4(Ca.x,Ca.y,Ca.z,0.0f);
   g_Cb[tid] = make_float4(Cb.x,Cb.y,Cb.z,0.0f);

   if (Vol_MType_K.y == 1)	// LE
   {
      tledSolverGPU_ROM_kernels::SPK_LE(tid,&Ca,&Cb,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 2)	// NH
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NH(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 3)	// TI
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TI(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 4)	// NHV
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NHV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 5)	// TIV
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TIV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 6)	// AB
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_AB(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 7)	// PY
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_PY(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
        
   // Save the stresses
   g_SPKa[tid] = make_float4(SPKa.x,SPKa.y,SPKa.z,0.0f);
   g_SPKb[tid] = make_float4(SPKb.x,SPKb.y,SPKb.z,0.0f);

   // Compute nodal forces: FEl = V*X*S*Dh'
   // Step 1: S = S*V
   SPKa = SPKa*Vol_MType_K.x;
   SPKb = SPKb*Vol_MType_K.x;
   // Step 2: A = X*S
   float A00 = XT00*SPKa.x+XT10*SPKb.x+XT20*SPKb.z;
   float A01 = XT00*SPKb.x+XT10*SPKa.y+XT20*SPKb.y;
   float A02 = XT00*SPKb.z+XT10*SPKb.y+XT20*SPKa.z;
   float A10 = XT01*SPKa.x+XT11*SPKb.x+XT21*SPKb.z;
   float A11 = XT01*SPKb.x+XT11*SPKa.y+XT21*SPKb.y;
   float A12 = XT01*SPKb.z+XT11*SPKb.y+XT21*SPKa.z;
   float A20 = XT02*SPKa.x+XT12*SPKb.x+XT22*SPKb.z;
   float A21 = XT02*SPKb.x+XT12*SPKa.y+XT22*SPKb.y;
   float A22 = XT02*SPKb.z+XT12*SPKb.y+XT22*SPKa.z;
   // Step 3: FEl = A*Dh'
   for (int i = 0; i < 4; i++)
   {
      float Dhi0 = Dh[i][0]; float Dhi1 = Dh[i][1]; float Dhi2 = Dh[i][2];
      tempF4.x = A00*Dhi0+A01*Dhi1+A02*Dhi2;
      tempF4.y = A10*Dhi0+A11*Dhi1+A12*Dhi2;
      tempF4.z = A20*Dhi0+A21*Dhi1+A22*Dhi2;

      g_FEl[tid*4 + i] = tempF4;
   }
} // tid < c_NumEls
}

// Compute internal forces for the current configuration - linear hexahedron
__global__ void ComputeNewForcesH8_kernel(float4* g_FEl, float4* g_StressStateIso, float4* g_StressStateVol,
                                          float4* g_SPKa, float4* g_SPKb, float4* g_Ca, float4* g_Cb)
{
   const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
if (tid < c_NumEls)
{
   float4 Vol_MType_K = tex1Dfetch(txVol_MType_K,tid);
   // Shape function derivatives
   float Dh[8][3];
   float4 tempF4;
   tempF4 = tex1Dfetch(txDhDx,6*tid);
   Dh[0][0] = tempF4.x; Dh[1][0] = tempF4.y; Dh[2][0] = tempF4.z; Dh[3][0] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,6*tid + 1);
   Dh[0][1] = tempF4.x; Dh[1][1] = tempF4.y; Dh[2][1] = tempF4.z; Dh[3][1] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,6*tid + 2);
   Dh[0][2] = tempF4.x; Dh[1][2] = tempF4.y; Dh[2][2] = tempF4.z; Dh[3][2] = tempF4.w;

   tempF4 = tex1Dfetch(txDhDx,6*tid + 3);
   Dh[4][0] = tempF4.x; Dh[5][0] = tempF4.y; Dh[6][0] = tempF4.z; Dh[7][0] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,6*tid + 4);
   Dh[4][1] = tempF4.x; Dh[5][1] = tempF4.y; Dh[6][1] = tempF4.z; Dh[7][1] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,6*tid + 5);
   Dh[4][2] = tempF4.x; Dh[5][2] = tempF4.y; Dh[6][2] = tempF4.z; Dh[7][2] = tempF4.w;

   // Displacements
   float u[8][3];
   int4 EInd;
   EInd = tex1Dfetch(txEInd,2*tid);
   tempF4 = tex1Dfetch(txUcurr,EInd.x);
   u[0][0] = tempF4.x; u[0][1] = tempF4.y; u[0][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.y);
   u[1][0] = tempF4.x; u[1][1] = tempF4.y; u[1][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.z);
   u[2][0] = tempF4.x; u[2][1] = tempF4.y; u[2][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.w);
   u[3][0] = tempF4.x; u[3][1] = tempF4.y; u[3][2] = tempF4.z;

   EInd = tex1Dfetch(txEInd,2*tid+1);
   tempF4 = tex1Dfetch(txUcurr,EInd.x);
   u[4][0] = tempF4.x; u[4][1] = tempF4.y; u[4][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.y);
   u[5][0] = tempF4.x; u[5][1] = tempF4.y; u[5][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.z);
   u[6][0] = tempF4.x; u[6][1] = tempF4.y; u[6][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.w);
   u[7][0] = tempF4.x; u[7][1] = tempF4.y; u[7][2] = tempF4.z;

   // // Deformation gradient (transpose of), X = Du + I, Du = u^T*Dh
   float XT00,XT01,XT02,XT10,XT11,XT12,XT20,XT21,XT22;
   XT00 = u[0][0]*Dh[0][0]+u[1][0]*Dh[1][0]+u[2][0]*Dh[2][0]+u[3][0]*Dh[3][0]+u[4][0]*Dh[4][0]+u[5][0]*Dh[5][0]+u[6][0]*Dh[6][0]+u[7][0]*Dh[7][0] + 1;
   XT01 = u[0][1]*Dh[0][0]+u[1][1]*Dh[1][0]+u[2][1]*Dh[2][0]+u[3][1]*Dh[3][0]+u[4][1]*Dh[4][0]+u[5][1]*Dh[5][0]+u[6][1]*Dh[6][0]+u[7][1]*Dh[7][0];
   XT02 = u[0][2]*Dh[0][0]+u[1][2]*Dh[1][0]+u[2][2]*Dh[2][0]+u[3][2]*Dh[3][0]+u[4][2]*Dh[4][0]+u[5][2]*Dh[5][0]+u[6][2]*Dh[6][0]+u[7][2]*Dh[7][0];
   XT10 = u[0][0]*Dh[0][1]+u[1][0]*Dh[1][1]+u[2][0]*Dh[2][1]+u[3][0]*Dh[3][1]+u[4][0]*Dh[4][1]+u[5][0]*Dh[5][1]+u[6][0]*Dh[6][1]+u[7][0]*Dh[7][1];
   XT11 = u[0][1]*Dh[0][1]+u[1][1]*Dh[1][1]+u[2][1]*Dh[2][1]+u[3][1]*Dh[3][1]+u[4][1]*Dh[4][1]+u[5][1]*Dh[5][1]+u[6][1]*Dh[6][1]+u[7][1]*Dh[7][1] + 1;
   XT12 = u[0][2]*Dh[0][1]+u[1][2]*Dh[1][1]+u[2][2]*Dh[2][1]+u[3][2]*Dh[3][1]+u[4][2]*Dh[4][1]+u[5][2]*Dh[5][1]+u[6][2]*Dh[6][1]+u[7][2]*Dh[7][1];
   XT20 = u[0][0]*Dh[0][2]+u[1][0]*Dh[1][2]+u[2][0]*Dh[2][2]+u[3][0]*Dh[3][2]+u[4][0]*Dh[4][2]+u[5][0]*Dh[5][2]+u[6][0]*Dh[6][2]+u[7][0]*Dh[7][2];
   XT21 = u[0][1]*Dh[0][2]+u[1][1]*Dh[1][2]+u[2][1]*Dh[2][2]+u[3][1]*Dh[3][2]+u[4][1]*Dh[4][2]+u[5][1]*Dh[5][2]+u[6][1]*Dh[6][2]+u[7][1]*Dh[7][2];
   XT22 = u[0][2]*Dh[0][2]+u[1][2]*Dh[1][2]+u[2][2]*Dh[2][2]+u[3][2]*Dh[3][2]+u[4][2]*Dh[4][2]+u[5][2]*Dh[5][2]+u[6][2]*Dh[6][2]+u[7][2]*Dh[7][2] + 1;

   // Stress
   float3 SPKa, SPKb;
   // Right Cauchy-Green deformation tensor
   float3 Ca = make_float3(XT00*XT00 + XT01*XT01 + XT02*XT02,	// C11
                     XT10*XT10 + XT11*XT11 + XT12*XT12,	// C22
                     XT20*XT20 + XT21*XT21 + XT22*XT22);	// C33
   float3 Cb = make_float3(XT00*XT10 + XT01*XT11 + XT02*XT12,	// C12
                     XT10*XT20 + XT11*XT21 + XT12*XT22,	// C23
                     XT00*XT20 + XT01*XT21 + XT02*XT22);	// C13
   
   // Save the deformation
   g_Ca[tid] = make_float4(Ca.x,Ca.y,Ca.z,0.0f);
   g_Cb[tid] = make_float4(Cb.x,Cb.y,Cb.z,0.0f);

   if (Vol_MType_K.y == 1)	// LE
   {
      tledSolverGPU_ROM_kernels::SPK_LE(tid,&Ca,&Cb,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 2)	// NH
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NH(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 3)	// TI
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TI(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 4)	// NHV
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NHV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 5)	// TIV
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TIV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 6)	// AB
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_AB(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 7)	// PY
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_PY(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   
   // Save the stresses
   g_SPKa[tid] = make_float4(SPKa.x,SPKa.y,SPKa.z,0.0f);
   g_SPKb[tid] = make_float4(SPKb.x,SPKb.y,SPKb.z,0.0f);

   // Compute nodal forces FEl = 8*detJ*BL'*S
   // (NB: For hex els this version seems to be faster than FEl = 8*detJ*X*S*Dh')
   SPKa = SPKa*Vol_MType_K.x;
   SPKb = SPKb*Vol_MType_K.x;
   float4 HG;
   for (int i = 0; i < 8; i++)
   {
      tempF4.x = Dh[i][0]*XT00*SPKa.x
               + Dh[i][1]*XT10*SPKa.y
               + Dh[i][2]*XT20*SPKa.z
               + (Dh[i][1]*XT00 + Dh[i][0]*XT10)*SPKb.x
               + (Dh[i][2]*XT10 + Dh[i][1]*XT20)*SPKb.y
               + (Dh[i][2]*XT00 + Dh[i][0]*XT20)*SPKb.z;
      tempF4.y = Dh[i][0]*XT01*SPKa.x
               + Dh[i][1]*XT11*SPKa.y
               + Dh[i][2]*XT21*SPKa.z
               + (Dh[i][1]*XT01 + Dh[i][0]*XT11)*SPKb.x
               + (Dh[i][2]*XT11 + Dh[i][1]*XT21)*SPKb.y
               + (Dh[i][2]*XT01 + Dh[i][0]*XT21)*SPKb.z;
      tempF4.z = Dh[i][0]*XT02*SPKa.x
               + Dh[i][1]*XT12*SPKa.y
               + Dh[i][2]*XT22*SPKa.z
               + (Dh[i][1]*XT02 + Dh[i][0]*XT12)*SPKb.x
               + (Dh[i][2]*XT12 + Dh[i][1]*XT22)*SPKb.y
               + (Dh[i][2]*XT02 + Dh[i][0]*XT22)*SPKb.z;

      // Hourglass control
      HG = tex1Dfetch(txHG,tid*16 + 2*i);
      tempF4.x += HG.x*u[0][0] + HG.y*u[1][0] + HG.z*u[2][0] + HG.w*u[3][0];	// Hourglass forces (= HG*u)
      tempF4.y += HG.x*u[0][1] + HG.y*u[1][1] + HG.z*u[2][1] + HG.w*u[3][1];
      tempF4.z += HG.x*u[0][2] + HG.y*u[1][2] + HG.z*u[2][2] + HG.w*u[3][2];
      HG = tex1Dfetch(txHG,tid*16 + 2*i + 1);
      tempF4.x += HG.x*u[4][0] + HG.y*u[5][0] + HG.z*u[6][0] + HG.w*u[7][0];
      tempF4.y += HG.x*u[4][1] + HG.y*u[5][1] + HG.z*u[6][1] + HG.w*u[7][1];
      tempF4.z += HG.x*u[4][2] + HG.y*u[5][2] + HG.z*u[6][2] + HG.w*u[7][2];

      g_FEl[tid*8 + i] = tempF4;
   }
} // (tid < g_NumEls)
}

// Compute effective loads for the current configuration
__global__ void ComputeEffectiveLoads_kernel(float4* g_Fint, float* g_f, tledGPUContacts* g_Contacts)
{
   const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
if (tid < c_NumNodes)
{
   // Gather internal forces
   float4 tempF4 = make_float4(0.0f,0.0f,0.0f,0.0f);
   int2 NodeMap = tex1Dfetch(txNodeMap,tid);
   int2 FCds;
   for (int i = 0; i < NodeMap.y; i++)
   {
      FCds = tex1Dfetch(txFCds,NodeMap.x + i);
      tempF4 = tempF4 + tex1Dfetch(txFEl,c_NPE*FCds.x + FCds.y);
   }
   g_Fint[tid] = tempF4;
   
   // Compute effective loads   
   tempF4 = tex1Dfetch(txF4Fext,tid) - tempF4; // Feff = Fext - Fint
   
   // Convert disp loads and contact disp into force loads, if applicable.
   // NB: Forces and displacements should not both be imposed on a given dof, nor
   // should contacted nodes be otherwise loaded
   // --> if this does happen, priority is 1. contacts, 2. disp loads, 3. forces
   int4 Mask = tex1Dfetch(txDispMask,tid);

   if (c_numContactCyls > 0)
   {
      float4 Ucurr = tex1Dfetch(txUcurr,tid);
      float4 Uprev = tex1Dfetch(txUprev,tid);
      float M = tex1Dfetch(txMfc,tid);
      
      float3 Unext = make_float3(Ucurr.x,Ucurr.y,Ucurr.z);
      float4 p = tex1Dfetch(txNodeCds,tid) + make_float4(Unext.x,Unext.y,Unext.z,0.0f);
      for (int i = 0; i < c_numContactCyls; i++)
      {
         int isSlv = g_Contacts->Cyls[i]->SlaveMask[tid];
         if (isSlv > 0) // Check if current node is a slave node
         {
            // NB: cylinder params arranged as
            // orig = (cyl_orig.x,cyl_orig.y,cyl_orig.z)
            // radius = cyl_orig.w
            // axis = (cyl_axis.x,cyl_axis.y,cyl_axis.z)
            // length = cyl_axis.w
   
            // Check penetration
            float4 cyl_orig = g_Contacts->Cyls[i]->Origin;
            float4 cyl_axis = g_Contacts->Cyls[i]->Axis;
            tledSolverGPU_ROM_kernels::ResolveCylinderContact(&p,&cyl_orig,&cyl_axis,&Unext);
         }
      }
      if (Unext.x != Ucurr.x)
         g_f[3*tid] = (Unext.x - c_gamma.y*Ucurr.x - c_gamma.z*Uprev.x)*M/c_gamma.x;
      else
         g_f[3*tid] = tempF4.x;
      if (Unext.y != Ucurr.y)
         g_f[3*tid+1] = (Unext.y - c_gamma.y*Ucurr.y - c_gamma.z*Uprev.y)*M/c_gamma.x;
      else
         g_f[3*tid+1] = tempF4.y;
      if (Unext.z != Ucurr.z)
         g_f[3*tid+2] = (Unext.z - c_gamma.y*Ucurr.z - c_gamma.z*Uprev.z)*M/c_gamma.x;
      else
         g_f[3*tid+2] = tempF4.z;
   }
   else if ( (Mask.x == 1)||(Mask.y == 1)||(Mask.z == 1) )
   {
      float4 Uload = tex1Dfetch(txUload,tid);
      float4 Ucurr = tex1Dfetch(txUcurr,tid);
      float4 Uprev = tex1Dfetch(txUprev,tid);
      float M = tex1Dfetch(txMfc,tid);
      
      if (Mask.x == 1)
         g_f[3*tid] = (Uload.x - c_gamma.y*Ucurr.x - c_gamma.z*Uprev.x)*M/c_gamma.x;
      else
         g_f[3*tid] = tempF4.x;
      if (Mask.y == 1)
         g_f[3*tid+1] = (Uload.y - c_gamma.y*Ucurr.y - c_gamma.z*Uprev.y)*M/c_gamma.x;
      else
         g_f[3*tid+1] = tempF4.y;
      if (Mask.z == 1)
         g_f[3*tid+2] = (Uload.z - c_gamma.y*Ucurr.z - c_gamma.z*Uprev.z)*M/c_gamma.x;
      else
         g_f[3*tid+2] = tempF4.z;
   }
   else
   {
      g_f[3*tid] = tempF4.x;
      g_f[3*tid+1] = tempF4.y;
      g_f[3*tid+2] = tempF4.z;
   }
      
} // tid < NumNodes
}

// Compute nodal displacements for the next step - central difference method
__global__ void ComputeNewDisps_kernel(float* g_UnextV, float4* g_Unext, bool* g_Divergence)
{
   const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
if (tid < c_NumNodes)
{
   float4 tempF4;
   // Unext acceleration term
   tempF4.x = g_UnextV[3*tid];
   tempF4.y = g_UnextV[3*tid+1];
   tempF4.z = g_UnextV[3*tid+2];
   float4 Unext = (c_gamma.x)*tempF4;
   // Add Ucurr term
   tempF4 = tex1Dfetch(txUcurr,tid);
   Unext = Unext + (c_gamma.y)*tempF4;
   // Add Uprev term
   tempF4 = tex1Dfetch(txUprev,tid);
   Unext = Unext + (c_gamma.z)*tempF4;
   // Get constraints
   int4 Mask = tex1Dfetch(txFixMask,tid);
   Unext = Unext*Mask;

   // FixMask (FM) is a binary mask indicating fixed DOFs
   // = 0 for fixed
   // = 1 otherwise

   // Check for divergence
   if (isnan(Unext.x+Unext.y+Unext.z))	// if any component is nan, the sum will be too
      *g_Divergence = true;

   g_Unext[tid] = Unext;

} // tid < NumNodes
}

// Compute element strain energy
__global__ void ComputeElementStrainEnergy_kernel(float* g_ElStrainEnergy, float4* g_SPKa, float4* g_SPKb, float4* g_Ca, float4* g_Cb)
{
   const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
if (tid < c_NumEls)
{
   // Element strain energy = (1/2) * int_V {Dot(Ev,Sv)} dV = (1/2) * V * Dot(Ev,Sv)
   
   float4 Vol_MType_K = tex1Dfetch(txVol_MType_K,tid);
   // Get first 1/2 of stress and strain vectors
   // Stresses
   float4 Sv = g_SPKa[tid];
   // Strains, E = (C-I)/2, Ev = [E00 E11 E22 2*E01 2*E12 2*E02]
   float4 Ev = g_Ca[tid];
   Ev.x = (Ev.x - 1)/2.0f; Ev.y = (Ev.y - 1)/2.0f; Ev.z = (Ev.z - 1)/2.0f;
   float e = Ev.x*Sv.x+Ev.y*Sv.y+Ev.z*Sv.z;
   // Get second half of stress and strain vectors
   Sv = g_SPKb[tid];
   Ev = g_Cb[tid];
   e += Ev.x*Sv.x+Ev.y*Sv.y+Ev.z*Sv.z;

   g_ElStrainEnergy[tid] = 0.5f*Vol_MType_K.x*e;
} // tid < NumEls
}

// Compute element stress for the current configuration - linear tetrahedron
__global__ void ComputeElementStressT4_kernel(float4* g_StressStateIso, float4* g_StressStateVol,
                                              float4* g_SPKa, float4* g_SPKb, float4* g_Ca, float4* g_Cb)
{
   const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
if (tid < c_NumEls)
{
   float4 Vol_MType_K = tex1Dfetch(txVol_MType_K,tid);
   // Shape function derivatives
   float Dh[4][3];
   float4 tempF4;
   tempF4 = tex1Dfetch(txDhDx,3*tid);
   Dh[0][0] = tempF4.x; Dh[1][0] = tempF4.y; Dh[2][0] = tempF4.z; Dh[3][0] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,3*tid + 1);
   Dh[0][1] = tempF4.x; Dh[1][1] = tempF4.y; Dh[2][1] = tempF4.z; Dh[3][1] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,3*tid + 2);
   Dh[0][2] = tempF4.x; Dh[1][2] = tempF4.y; Dh[2][2] = tempF4.z; Dh[3][2] = tempF4.w;

   // Displacements
   float u[4][3];
   int4 EInd;
   EInd = tex1Dfetch(txEInd,tid);
   tempF4 = tex1Dfetch(txUcurr,EInd.x);
   u[0][0] = tempF4.x; u[0][1] = tempF4.y; u[0][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.y);
   u[1][0] = tempF4.x; u[1][1] = tempF4.y; u[1][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.z);
   u[2][0] = tempF4.x; u[2][1] = tempF4.y; u[2][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.w);
   u[3][0] = tempF4.x; u[3][1] = tempF4.y; u[3][2] = tempF4.z;

   // Deformation gradient (transpose of), X = Du + I, Du = u^T*Dh
   float XT00,XT01,XT02,XT10,XT11,XT12,XT20,XT21,XT22;
   XT00 = u[0][0]*Dh[0][0]+u[1][0]*Dh[1][0]+u[2][0]*Dh[2][0]+u[3][0]*Dh[3][0] + 1;
   XT01 = u[0][1]*Dh[0][0]+u[1][1]*Dh[1][0]+u[2][1]*Dh[2][0]+u[3][1]*Dh[3][0];
   XT02 = u[0][2]*Dh[0][0]+u[1][2]*Dh[1][0]+u[2][2]*Dh[2][0]+u[3][2]*Dh[3][0];
   XT10 = u[0][0]*Dh[0][1]+u[1][0]*Dh[1][1]+u[2][0]*Dh[2][1]+u[3][0]*Dh[3][1];
   XT11 = u[0][1]*Dh[0][1]+u[1][1]*Dh[1][1]+u[2][1]*Dh[2][1]+u[3][1]*Dh[3][1] + 1;
   XT12 = u[0][2]*Dh[0][1]+u[1][2]*Dh[1][1]+u[2][2]*Dh[2][1]+u[3][2]*Dh[3][1];
   XT20 = u[0][0]*Dh[0][2]+u[1][0]*Dh[1][2]+u[2][0]*Dh[2][2]+u[3][0]*Dh[3][2];
   XT21 = u[0][1]*Dh[0][2]+u[1][1]*Dh[1][2]+u[2][1]*Dh[2][2]+u[3][1]*Dh[3][2];
   XT22 = u[0][2]*Dh[0][2]+u[1][2]*Dh[1][2]+u[2][2]*Dh[2][2]+u[3][2]*Dh[3][2] + 1;

   // Stress
   float3 SPKa, SPKb;
   // Right Cauchy-Green deformation tensor
   float3 Ca = make_float3(XT00*XT00 + XT01*XT01 + XT02*XT02,   // C11
                           XT10*XT10 + XT11*XT11 + XT12*XT12,   // C22
                           XT20*XT20 + XT21*XT21 + XT22*XT22);  // C33
   float3 Cb = make_float3(XT00*XT10 + XT01*XT11 + XT02*XT12,   // C12
                           XT10*XT20 + XT11*XT21 + XT12*XT22,   // C23
                           XT00*XT20 + XT01*XT21 + XT02*XT22);  // C13
        
   // Save the deformation
   g_Ca[tid] = make_float4(Ca.x,Ca.y,Ca.z,0.0f);
   g_Cb[tid] = make_float4(Cb.x,Cb.y,Cb.z,0.0f);

   if (Vol_MType_K.y == 1)	// LE
   {
      tledSolverGPU_ROM_kernels::SPK_LE(tid,&Ca,&Cb,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 2)	// NH
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NH(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 3)	// TI
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TI(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 4)	// NHV
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NHV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 5)	// TIV
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TIV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 6)	// AB
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_AB(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 7)	// PY
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_PY(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
        
   // Save the stresses
   g_SPKa[tid] = make_float4(SPKa.x,SPKa.y,SPKa.z,0.0f);
   g_SPKb[tid] = make_float4(SPKb.x,SPKb.y,SPKb.z,0.0f);
} // tid < NumEls
}

// Compute element stress for the current configuration - linear hexahedron
__global__ void ComputeElementStressH8_kernel(float4* g_StressStateIso, float4* g_StressStateVol,
                                              float4* g_SPKa, float4* g_SPKb, float4* g_Ca, float4* g_Cb)
{
   const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
if (tid < c_NumEls)
{
   float4 Vol_MType_K = tex1Dfetch(txVol_MType_K,tid);
   // Shape function derivatives
   float Dh[8][3];
   float4 tempF4;
   tempF4 = tex1Dfetch(txDhDx,6*tid);
   Dh[0][0] = tempF4.x; Dh[1][0] = tempF4.y; Dh[2][0] = tempF4.z; Dh[3][0] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,6*tid + 1);
   Dh[0][1] = tempF4.x; Dh[1][1] = tempF4.y; Dh[2][1] = tempF4.z; Dh[3][1] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,6*tid + 2);
   Dh[0][2] = tempF4.x; Dh[1][2] = tempF4.y; Dh[2][2] = tempF4.z; Dh[3][2] = tempF4.w;

   tempF4 = tex1Dfetch(txDhDx,6*tid + 3);
   Dh[4][0] = tempF4.x; Dh[5][0] = tempF4.y; Dh[6][0] = tempF4.z; Dh[7][0] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,6*tid + 4);
   Dh[4][1] = tempF4.x; Dh[5][1] = tempF4.y; Dh[6][1] = tempF4.z; Dh[7][1] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,6*tid + 5);
   Dh[4][2] = tempF4.x; Dh[5][2] = tempF4.y; Dh[6][2] = tempF4.z; Dh[7][2] = tempF4.w;

   // Displacements
   float u[8][3];
   int4 EInd;
   EInd = tex1Dfetch(txEInd,2*tid);
   tempF4 = tex1Dfetch(txUcurr,EInd.x);
   u[0][0] = tempF4.x; u[0][1] = tempF4.y; u[0][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.y);
   u[1][0] = tempF4.x; u[1][1] = tempF4.y; u[1][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.z);
   u[2][0] = tempF4.x; u[2][1] = tempF4.y; u[2][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.w);
   u[3][0] = tempF4.x; u[3][1] = tempF4.y; u[3][2] = tempF4.z;

   EInd = tex1Dfetch(txEInd,2*tid+1);
   tempF4 = tex1Dfetch(txUcurr,EInd.x);
   u[4][0] = tempF4.x; u[4][1] = tempF4.y; u[4][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.y);
   u[5][0] = tempF4.x; u[5][1] = tempF4.y; u[5][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.z);
   u[6][0] = tempF4.x; u[6][1] = tempF4.y; u[6][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.w);
   u[7][0] = tempF4.x; u[7][1] = tempF4.y; u[7][2] = tempF4.z;

   // // Deformation gradient (transpose of), X = Du + I, Du = u^T*Dh
   float XT00,XT01,XT02,XT10,XT11,XT12,XT20,XT21,XT22;
   XT00 = u[0][0]*Dh[0][0]+u[1][0]*Dh[1][0]+u[2][0]*Dh[2][0]+u[3][0]*Dh[3][0]+u[4][0]*Dh[4][0]+u[5][0]*Dh[5][0]+u[6][0]*Dh[6][0]+u[7][0]*Dh[7][0] + 1;
   XT01 = u[0][1]*Dh[0][0]+u[1][1]*Dh[1][0]+u[2][1]*Dh[2][0]+u[3][1]*Dh[3][0]+u[4][1]*Dh[4][0]+u[5][1]*Dh[5][0]+u[6][1]*Dh[6][0]+u[7][1]*Dh[7][0];
   XT02 = u[0][2]*Dh[0][0]+u[1][2]*Dh[1][0]+u[2][2]*Dh[2][0]+u[3][2]*Dh[3][0]+u[4][2]*Dh[4][0]+u[5][2]*Dh[5][0]+u[6][2]*Dh[6][0]+u[7][2]*Dh[7][0];
   XT10 = u[0][0]*Dh[0][1]+u[1][0]*Dh[1][1]+u[2][0]*Dh[2][1]+u[3][0]*Dh[3][1]+u[4][0]*Dh[4][1]+u[5][0]*Dh[5][1]+u[6][0]*Dh[6][1]+u[7][0]*Dh[7][1];
   XT11 = u[0][1]*Dh[0][1]+u[1][1]*Dh[1][1]+u[2][1]*Dh[2][1]+u[3][1]*Dh[3][1]+u[4][1]*Dh[4][1]+u[5][1]*Dh[5][1]+u[6][1]*Dh[6][1]+u[7][1]*Dh[7][1] + 1;
   XT12 = u[0][2]*Dh[0][1]+u[1][2]*Dh[1][1]+u[2][2]*Dh[2][1]+u[3][2]*Dh[3][1]+u[4][2]*Dh[4][1]+u[5][2]*Dh[5][1]+u[6][2]*Dh[6][1]+u[7][2]*Dh[7][1];
   XT20 = u[0][0]*Dh[0][2]+u[1][0]*Dh[1][2]+u[2][0]*Dh[2][2]+u[3][0]*Dh[3][2]+u[4][0]*Dh[4][2]+u[5][0]*Dh[5][2]+u[6][0]*Dh[6][2]+u[7][0]*Dh[7][2];
   XT21 = u[0][1]*Dh[0][2]+u[1][1]*Dh[1][2]+u[2][1]*Dh[2][2]+u[3][1]*Dh[3][2]+u[4][1]*Dh[4][2]+u[5][1]*Dh[5][2]+u[6][1]*Dh[6][2]+u[7][1]*Dh[7][2];
   XT22 = u[0][2]*Dh[0][2]+u[1][2]*Dh[1][2]+u[2][2]*Dh[2][2]+u[3][2]*Dh[3][2]+u[4][2]*Dh[4][2]+u[5][2]*Dh[5][2]+u[6][2]*Dh[6][2]+u[7][2]*Dh[7][2] + 1;

   // Stress
   float3 SPKa, SPKb;
   // Right Cauchy-Green deformation tensor
   float3 Ca = make_float3(XT00*XT00 + XT01*XT01 + XT02*XT02,	// C11
                           XT10*XT10 + XT11*XT11 + XT12*XT12,	// C22
                           XT20*XT20 + XT21*XT21 + XT22*XT22);	// C33
   float3 Cb = make_float3(XT00*XT10 + XT01*XT11 + XT02*XT12,	// C12
                           XT10*XT20 + XT11*XT21 + XT12*XT22,	// C23
                           XT00*XT20 + XT01*XT21 + XT02*XT22);	// C13
   
   // Save the deformation
   g_Ca[tid] = make_float4(Ca.x,Ca.y,Ca.z,0.0f);
   g_Cb[tid] = make_float4(Cb.x,Cb.y,Cb.z,0.0f);

   if (Vol_MType_K.y == 1)	// LE
   {
      tledSolverGPU_ROM_kernels::SPK_LE(tid,&Ca,&Cb,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 2)	// NH
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NH(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 3)	// TI
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TI(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 4)	// NHV
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NHV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 5)	// TIV
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TIV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 6)	// AB
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_AB(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 7)	// PY
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_PY(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   
   // Save the stresses
   g_SPKa[tid] = make_float4(SPKa.x,SPKa.y,SPKa.z,0.0f);
   g_SPKb[tid] = make_float4(SPKb.x,SPKb.y,SPKb.z,0.0f);
} // tid < NumEls
}

// Compute element stress for the current configuration - nodal-averaged pressure tetrahedron
__global__ void ComputeModifiedElementStressT4ANP_kernel(float4* g_StressStateIso, float4* g_StressStateVol,
      float4* g_SPKa, float4* g_SPKb, float4* g_Ca, float4* g_Cb)
{
   const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
if (tid < c_NumEls)
{
   float4 Vol_MType_K = tex1Dfetch(txVol_MType_K,tid);	// Vol = Vol_MType_K.x, MType = Vol_MType_K.y, K = Vol_MType_K.z
   // Shape function derivatives
   float Dh[4][3];
   float4 tempF4;
   tempF4 = tex1Dfetch(txDhDx,3*tid);
   Dh[0][0] = tempF4.x; Dh[1][0] = tempF4.y; Dh[2][0] = tempF4.z; Dh[3][0] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,3*tid + 1);
   Dh[0][1] = tempF4.x; Dh[1][1] = tempF4.y; Dh[2][1] = tempF4.z; Dh[3][1] = tempF4.w;
   tempF4 = tex1Dfetch(txDhDx,3*tid + 2);
   Dh[0][2] = tempF4.x; Dh[1][2] = tempF4.y; Dh[2][2] = tempF4.z; Dh[3][2] = tempF4.w;

   // Displacements
   float u[4][3];
   int4 EInd;
   EInd = tex1Dfetch(txEInd,tid);
   tempF4 = tex1Dfetch(txUcurr,EInd.x);
   u[0][0] = tempF4.x; u[0][1] = tempF4.y; u[0][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.y);
   u[1][0] = tempF4.x; u[1][1] = tempF4.y; u[1][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.z);
   u[2][0] = tempF4.x; u[2][1] = tempF4.y; u[2][2] = tempF4.z;
   tempF4 = tex1Dfetch(txUcurr,EInd.w);
   u[3][0] = tempF4.x; u[3][1] = tempF4.y; u[3][2] = tempF4.z;

   // Deformation gradient (transpose of), X = Du + I, Du = u^T*Dh
   float XT00,XT01,XT02,XT10,XT11,XT12,XT20,XT21,XT22;
   XT00 = u[0][0]*Dh[0][0]+u[1][0]*Dh[1][0]+u[2][0]*Dh[2][0]+u[3][0]*Dh[3][0] + 1;
   XT01 = u[0][1]*Dh[0][0]+u[1][1]*Dh[1][0]+u[2][1]*Dh[2][0]+u[3][1]*Dh[3][0];
   XT02 = u[0][2]*Dh[0][0]+u[1][2]*Dh[1][0]+u[2][2]*Dh[2][0]+u[3][2]*Dh[3][0];
   XT10 = u[0][0]*Dh[0][1]+u[1][0]*Dh[1][1]+u[2][0]*Dh[2][1]+u[3][0]*Dh[3][1];
   XT11 = u[0][1]*Dh[0][1]+u[1][1]*Dh[1][1]+u[2][1]*Dh[2][1]+u[3][1]*Dh[3][1] + 1;
   XT12 = u[0][2]*Dh[0][1]+u[1][2]*Dh[1][1]+u[2][2]*Dh[2][1]+u[3][2]*Dh[3][1];
   XT20 = u[0][0]*Dh[0][2]+u[1][0]*Dh[1][2]+u[2][0]*Dh[2][2]+u[3][0]*Dh[3][2];
   XT21 = u[0][1]*Dh[0][2]+u[1][1]*Dh[1][2]+u[2][1]*Dh[2][2]+u[3][1]*Dh[3][2];
   XT22 = u[0][2]*Dh[0][2]+u[1][2]*Dh[1][2]+u[2][2]*Dh[2][2]+u[3][2]*Dh[3][2] + 1;

   float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
   J = powf(J,-1.0f/3.0f); // J = J^-1/3

   // Compute modified pressure
   float Pm = 0.0f;
   Pm += tex1Dfetch(txPa,EInd.x);
   Pm += tex1Dfetch(txPa,EInd.y);
   Pm += tex1Dfetch(txPa,EInd.z);
   Pm += tex1Dfetch(txPa,EInd.w);
   Pm /= 16;

   // Compute modified Jacobian
   float Jm = Pm/Vol_MType_K.z + 1; // Jm = Pm/K + 1
   // Compute modified def. grad.
   Jm = powf(Jm,1.0f/3.0f); // Jm = Jm^1/3
   Jm = Jm*J;	// Jm = Jm^1/3 * J^-1/3
   XT00 *= Jm;
   XT01 *= Jm;
   XT02 *= Jm;
   XT10 *= Jm;
   XT11 *= Jm;
   XT12 *= Jm;
   XT20 *= Jm;
   XT21 *= Jm;
   XT22 *= Jm;

   // Stress
   float3 SPKa, SPKb;
   // Right Cauchy-Green deformation tensor
   float3 Ca = make_float3(XT00*XT00 + XT01*XT01 + XT02*XT02,	// C11
                           XT10*XT10 + XT11*XT11 + XT12*XT12,	// C22
                           XT20*XT20 + XT21*XT21 + XT22*XT22);	// C33
   float3 Cb = make_float3(XT00*XT10 + XT01*XT11 + XT02*XT12,	// C12
                           XT10*XT20 + XT11*XT21 + XT12*XT22,	// C23
                           XT00*XT20 + XT01*XT21 + XT02*XT22);	// C13
        
   // Save the deformation
   g_Ca[tid] = make_float4(Ca.x,Ca.y,Ca.z,0.0f);
   g_Cb[tid] = make_float4(Cb.x,Cb.y,Cb.z,0.0f);

   if (Vol_MType_K.y == 1)	// LE
   {
      tledSolverGPU_ROM_kernels::SPK_LE(tid,&Ca,&Cb,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 2)	// NH
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NH(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 3)	// TI
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TI(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 4)	// NHV
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_NHV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 5)	// TIV
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_TIV(tid,&Ca,&Cb,J,g_StressStateIso,g_StressStateVol,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 6)	// AB
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_AB(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
   else if (Vol_MType_K.y == 7)	// PY
   {
      float J = XT00*(XT11*XT22 - XT21*XT12) + XT01*(XT12*XT20 - XT10*XT22) + XT02*(XT10*XT21 - XT11*XT20);
      tledSolverGPU_ROM_kernels::SPK_PY(tid,&Ca,&Cb,J,&SPKa,&SPKb);
   }
        
   // Save the stresses
   g_SPKa[tid] = make_float4(SPKa.x,SPKa.y,SPKa.z,0.0f);
   g_SPKb[tid] = make_float4(SPKb.x,SPKb.y,SPKb.z,0.0f);
} // tid < c_NumEls
}


} // namespace tledSolverGPU_ROM_kernels



#endif // tledSolverGPU_ROM_kernels_CU
