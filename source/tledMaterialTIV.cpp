// =========================================================================
// File:       tledMaterialTIV.cpp
// Purpose:    Material class - Transversely isotropic neo-Hookean viscoelastic
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    January 2008
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#include "tledMaterialTIV.h"
#include "tledMatrixFunctions.h"

#include <iostream>
#include <cmath>

using namespace std;

tledMaterialTIV::tledMaterialTIV()
{
}

tledMaterialTIV::tledMaterialTIV(float* MatParams, float Density) : tledMaterial(Density)
{   
   M = MatParams[0];
   K = MatParams[1];
   E = MatParams[2];
   float a[3] = {MatParams[3],MatParams[4],MatParams[5]};
   float norma = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
   if (norma != 1.0f)
   {
      // Normalise the direction vector
      for (int i = 0; i < 3; i++)
         a[i] /= norma;
   }
   A11 = a[0]*a[0];	// A(1,1) = a1*a1
   A22 = a[1]*a[1];	// A(2,2) = a2*a2
   A33 = a[2]*a[2];	// A(3,3) = a3*a3
   A12 = a[0]*a[1];	// A(1,2) = a1*a2
   A23 = a[1]*a[2];	// A(2,3) = a2*a3
   A13 = a[0]*a[2];	// A(1,3) = a1*a3
   Dt = MatParams[6];
   Ni = (int)MatParams[7];
   Nv = (int)MatParams[8];
   // Set up isochoric terms
   Ai = new float[2*Ni];
   for (int i = 0; i < Ni; i++)
   {
      Ai[2*i]   = Dt*MatParams[2*i+9]/(Dt + MatParams[2*i+10]);
      Ai[2*i+1] = MatParams[2*i+10]/(Dt + MatParams[2*i+10]);
   }
   Di = new float[6*Ni];
   memset(Di,0,sizeof(float)*6*Ni);
   // Set up volumetric terms
   Av = new float[2*Nv];
   for (int i = 0; i < Nv; i++)
   {
      Av[2*i]   = Dt*MatParams[2*i+9+2*Ni]/(Dt + MatParams[2*i+10+2*Ni]);
      Av[2*i+1] = MatParams[2*i+10+2*Ni]/(Dt + MatParams[2*i+10+2*Ni]);
   }
   Dv = new float[6*Nv];
   memset(Dv,0,sizeof(float)*6*Nv);
}

tledMaterialTIV::~tledMaterialTIV()
{
   delete Ai;
   delete Av;
   delete Di;
   delete Dv;
}

void tledMaterialTIV::GetHGLameParams(float* Lambda, float* Mu)
{
   *Lambda = K - 2*M/3;
   *Mu = M;
}

void tledMaterialTIV::SetTimeStep(float dt)
{
   float Dt_old = Dt;
   Dt = dt;
   // Update visc integration coeffs
   for (int i = 0; i < Ni; i++) // Iso terms
   {
      float A = Ai[2*i];
      float B = Ai[2*i+1];
      float t = Dt_old/(1/B-1);
      float a = A*(Dt_old+t)/Dt_old;
      Ai[2*i] = Dt*a/(Dt+t); // New vals
      Ai[2*i+1] = t/(Dt+t);
   }
   for (int i = 0; i < Nv; i++) // Vol terms
   {
      float A = Av[2*i];
      float B = Av[2*i+1];
      float t = Dt_old/(1/B-1);
      float a = A*(Dt_old+t)/Dt_old;
      Av[2*i] = Dt*a/(Dt+t); // New vals
      Av[2*i+1] = t/(Dt+t);
   }
}

void tledMaterialTIV::SetParams(vector<float> params)
{
   M = params[0];
   K = params[1];
   E = params[2];
   float a[3] = {params[3],params[4],params[5]};
   float norma = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
   if (norma != 1.0f)
   {
      // Normalise the direction vector
      for (int i = 0; i < 3; i++)
         a[i] /= norma;
   }
   A11 = a[0]*a[0];	// A(1,1) = a1*a1
   A22 = a[1]*a[1];	// A(2,2) = a2*a2
   A33 = a[2]*a[2];	// A(3,3) = a3*a3
   A12 = a[0]*a[1];	// A(1,2) = a1*a2
   A23 = a[1]*a[2];	// A(2,3) = a2*a3
   A13 = a[0]*a[2];	// A(1,3) = a1*a3
   Dt = params[6];
   Ni = (int)params[7];
   Nv = (int)params[8];
   // Set up isochoric terms
   delete Ai;
   Ai = new float[2*Ni];
   for (int i = 0; i < Ni; i++)
   {
      Ai[2*i]   = Dt*params[2*i+9]/(Dt + params[2*i+10]);
      Ai[2*i+1] = params[2*i+10]/(Dt + params[2*i+10]);
   }
   delete Di;
   Di = new float[6*Ni];
   memset(Di,0,sizeof(float)*6*Ni);
   // Set up volumetric terms
   delete Av;
   Av = new float[2*Nv];
   for (int i = 0; i < Nv; i++)
   {
      Av[2*i]   = Dt*params[2*i+9+2*Ni]/(Dt + params[2*i+10+2*Ni]);
      Av[2*i+1] = params[2*i+10+2*Ni]/(Dt + params[2*i+10+2*Ni]);
   }
   delete Dv;
   Dv = new float[6*Nv];
   memset(Dv,0,sizeof(float)*6*Nv);
}

void tledMaterialTIV::ComputeSPKStress(float X[3][3], float *SPK)
{
   int i;
   float J,J23,invdetC,x1,x2,x3,x4,x5;
   float XT11,XT12,XT13,XT21,XT22,XT23,XT31,XT32,XT33;
   float C11,C12,C13,C22,C23,C33;
   float Ci11,Ci12,Ci13,Ci22,Ci23,Ci33;
   float SiE11,SiE12,SiE13,SiE22,SiE23,SiE33;
   float SvE11,SvE12,SvE13,SvE22,SvE23,SvE33;

   // Transpose of deformation gradient
   XT11 = X[0][0]; XT12 = X[1][0]; XT13 = X[2][0];
   XT21 = X[0][1]; XT22 = X[1][1]; XT23 = X[2][1];
   XT31 = X[0][2]; XT32 = X[1][2]; XT33 = X[2][2];

   // Right Cauchy-Green deformation tensor
   C11 = XT11*XT11 + XT12*XT12 + XT13*XT13;
   C12 = XT11*XT21 + XT12*XT22 + XT13*XT23;
   C13 = XT11*XT31 + XT12*XT32 + XT13*XT33;
   C22 = XT21*XT21 + XT22*XT22 + XT23*XT23;
   C23 = XT21*XT31 + XT22*XT32 + XT23*XT33;
   C33 = XT31*XT31 + XT32*XT32 + XT33*XT33;

   // Determinant of deformation gradient
   J = XT11*(XT22*XT33 - XT32*XT23) + XT12*(XT23*XT31 - XT21*XT33) + XT13*(XT21*XT32 - XT22*XT31);
   J23 = (float)pow((double)J,-(double)2/(double)3);	// J23 = J^(-2/3)
   x1 = J23*M;
   x2 = J23*(A11*C11+A22*C22+A33*C33+2*A12*C12+2*A23*C23+2*A13*C13) - 1; // Bracketed term is I4 = A:C
   x3 = J23*E*x2;
   x4 = -(E*x2*(x2+1) + x1*(C11+C22+C33))/3;
   x5 = K*J*(J-1);

   // C inverse
   invdetC = 1/(C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23));
   Ci11 = (C22*C33 - C23*C23)*invdetC;
   Ci12 = (C13*C23 - C12*C33)*invdetC;
   Ci13 = (C12*C23 - C13*C22)*invdetC;
   Ci22 = (C11*C33 - C13*C13)*invdetC;
   Ci23 = (C12*C13 - C11*C23)*invdetC;
   Ci33 = (C11*C22 - C12*C12)*invdetC;

   // Elastic components of response *********************************
   // SiE
   SiE11 = x3*A11 + x4*Ci11 + x1;
   SiE22 = x3*A22 + x4*Ci22 + x1;
   SiE33 = x3*A33 + x4*Ci33 + x1;
   SiE12 = x3*A12 + x4*Ci12;
   SiE23 = x3*A23 + x4*Ci23;
   SiE13 = x3*A13 + x4*Ci13;
   // SvE
   SvE11 = x5*Ci11;
   SvE22 = x5*Ci22;
   SvE33 = x5*Ci33;
   SvE12 = x5*Ci12;
   SvE23 = x5*Ci23;
   SvE13 = x5*Ci13;

   // State variables ************************************************
   // Reuse C-variables for summing state variables
   // --> not required anymore, and saves declaring another 6 floats for this purpose
   C11 = 0; C22 = 0; C33 = 0; C12 = 0; C23 = 0; C13 = 0;
   // Isochoric
   for (i = 0; i < Ni; i++)
   {
      Di[i*6]   *= Ai[2*i+1]; Di[i*6]   += Ai[2*i]*SiE11; C11 += Di[i*6];
      Di[i*6+1] *= Ai[2*i+1]; Di[i*6+1] += Ai[2*i]*SiE22; C22 += Di[i*6+1];
      Di[i*6+2] *= Ai[2*i+1]; Di[i*6+2] += Ai[2*i]*SiE33; C33 += Di[i*6+2];
      Di[i*6+3] *= Ai[2*i+1]; Di[i*6+3] += Ai[2*i]*SiE12; C12 += Di[i*6+3];
      Di[i*6+4] *= Ai[2*i+1]; Di[i*6+4] += Ai[2*i]*SiE23; C23 += Di[i*6+4];
      Di[i*6+5] *= Ai[2*i+1]; Di[i*6+5] += Ai[2*i]*SiE13; C13 += Di[i*6+5];
   }
   // Volumetric
   for (i = 0; i < Nv; i++)
   {
      Dv[i*6]   *= Av[2*i+1]; Dv[i*6]   += Av[2*i]*SvE11; C11 += Dv[i*6];
      Dv[i*6+1] *= Av[2*i+1]; Dv[i*6+1] += Av[2*i]*SvE22; C22 += Dv[i*6+1];
      Dv[i*6+2] *= Av[2*i+1]; Dv[i*6+2] += Av[2*i]*SvE33; C33 += Dv[i*6+2];
      Dv[i*6+3] *= Av[2*i+1]; Dv[i*6+3] += Av[2*i]*SvE12; C12 += Dv[i*6+3];
      Dv[i*6+4] *= Av[2*i+1]; Dv[i*6+4] += Av[2*i]*SvE23; C23 += Dv[i*6+4];
      Dv[i*6+5] *= Av[2*i+1]; Dv[i*6+5] += Av[2*i]*SvE13; C13 += Dv[i*6+5];
   }

   // Total stress ***************************************************
   SPK[0] = SiE11 + SvE11 - C11;
   SPK[1] = SiE22 + SvE22 - C22;
   SPK[2] = SiE33 + SvE33 - C33;
   SPK[3] = SiE12 + SvE12 - C12;
   SPK[4] = SiE23 + SvE23 - C23;
   SPK[5] = SiE13 + SvE13 - C13;
}

void tledMaterialTIV::ComputeSPKStress(float X[3][3], float SPK[3][3])
{
   int i;
   float J,J23,invdetC,x1,x2,x3,x4,x5;
   float XT11,XT12,XT13,XT21,XT22,XT23,XT31,XT32,XT33;
   float C11,C12,C13,C22,C23,C33;
   float Ci11,Ci12,Ci13,Ci22,Ci23,Ci33;
   float SiE11,SiE12,SiE13,SiE22,SiE23,SiE33;
   float SvE11,SvE12,SvE13,SvE22,SvE23,SvE33;

   // Transpose of deformation gradient
   XT11 = X[0][0]; XT12 = X[1][0]; XT13 = X[2][0];
   XT21 = X[0][1]; XT22 = X[1][1]; XT23 = X[2][1];
   XT31 = X[0][2]; XT32 = X[1][2]; XT33 = X[2][2];

   // Right Cauchy-Green deformation tensor
   C11 = XT11*XT11 + XT12*XT12 + XT13*XT13;
   C12 = XT11*XT21 + XT12*XT22 + XT13*XT23;
   C13 = XT11*XT31 + XT12*XT32 + XT13*XT33;
   C22 = XT21*XT21 + XT22*XT22 + XT23*XT23;
   C23 = XT21*XT31 + XT22*XT32 + XT23*XT33;
   C33 = XT31*XT31 + XT32*XT32 + XT33*XT33;

   // Determinant of deformation gradient
   J = XT11*(XT22*XT33 - XT32*XT23) + XT12*(XT23*XT31 - XT21*XT33) + XT13*(XT21*XT32 - XT22*XT31);
   J23 = (float)pow((double)J,-(double)2/(double)3);	// J23 = J^(-2/3)
   x1 = J23*M;
   x2 = J23*(A11*C11+A22*C22+A33*C33+2*A12*C12+2*A23*C23+2*A13*C13) - 1; // Bracketed term is I4 = A:C
   x3 = J23*E*x2;
   x4 = -(E*x2*(x2+1) + x1*(C11+C22+C33))/3;
   x5 = K*J*(J-1);

   // C inverse
   invdetC = 1/(C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23));
   Ci11 = (C22*C33 - C23*C23)*invdetC;
   Ci12 = (C13*C23 - C12*C33)*invdetC;
   Ci13 = (C12*C23 - C13*C22)*invdetC;
   Ci22 = (C11*C33 - C13*C13)*invdetC;
   Ci23 = (C12*C13 - C11*C23)*invdetC;
   Ci33 = (C11*C22 - C12*C12)*invdetC;

   // Elastic components of response *********************************
   // SiE
   SiE11 = x3*A11 + x4*Ci11 + x1;
   SiE22 = x3*A22 + x4*Ci22 + x1;
   SiE33 = x3*A33 + x4*Ci33 + x1;
   SiE12 = x3*A12 + x4*Ci12;
   SiE23 = x3*A23 + x4*Ci23;
   SiE13 = x3*A13 + x4*Ci13;
   // SvE
   SvE11 = x5*Ci11;
   SvE22 = x5*Ci22;
   SvE33 = x5*Ci33;
   SvE12 = x5*Ci12;
   SvE23 = x5*Ci23;
   SvE13 = x5*Ci13;

   // State variables ************************************************
   // Reuse C-variables for summing state variables
   // --> not required anymore, and saves declaring another 6 floats for this purpose
   C11 = 0; C22 = 0; C33 = 0; C12 = 0; C23 = 0; C13 = 0;
   // Isochoric
   for (i = 0; i < Ni; i++)
   {
      Di[i*6]   *= Ai[2*i+1]; Di[i*6]   += Ai[2*i]*SiE11; C11 += Di[i*6];
      Di[i*6+1] *= Ai[2*i+1]; Di[i*6+1] += Ai[2*i]*SiE22; C22 += Di[i*6+1];
      Di[i*6+2] *= Ai[2*i+1]; Di[i*6+2] += Ai[2*i]*SiE33; C33 += Di[i*6+2];
      Di[i*6+3] *= Ai[2*i+1]; Di[i*6+3] += Ai[2*i]*SiE12; C12 += Di[i*6+3];
      Di[i*6+4] *= Ai[2*i+1]; Di[i*6+4] += Ai[2*i]*SiE23; C23 += Di[i*6+4];
      Di[i*6+5] *= Ai[2*i+1]; Di[i*6+5] += Ai[2*i]*SiE13; C13 += Di[i*6+5];
   }
   // Volumetric
   for (i = 0; i < Nv; i++)
   {
      Dv[i*6]   *= Av[2*i+1]; Dv[i*6]   += Av[2*i]*SvE11; C11 += Dv[i*6];
      Dv[i*6+1] *= Av[2*i+1]; Dv[i*6+1] += Av[2*i]*SvE22; C22 += Dv[i*6+1];
      Dv[i*6+2] *= Av[2*i+1]; Dv[i*6+2] += Av[2*i]*SvE33; C33 += Dv[i*6+2];
      Dv[i*6+3] *= Av[2*i+1]; Dv[i*6+3] += Av[2*i]*SvE12; C12 += Dv[i*6+3];
      Dv[i*6+4] *= Av[2*i+1]; Dv[i*6+4] += Av[2*i]*SvE23; C23 += Dv[i*6+4];
      Dv[i*6+5] *= Av[2*i+1]; Dv[i*6+5] += Av[2*i]*SvE13; C13 += Dv[i*6+5];
   }

   // Total stress ***************************************************
   SPK[0][0] = SiE11 + SvE11 - C11; // S00
   SPK[1][1] = SiE22 + SvE22 - C22; // S11
   SPK[2][2] = SiE33 + SvE33 - C33; // S22
   SPK[0][1] = SiE12 + SvE12 - C12; // S01
   SPK[1][0] = SPK[0][1];
   SPK[1][2] = SiE23 + SvE23 - C23; // S12
   SPK[2][1] = SPK[1][2];
   SPK[0][2] = SiE13 + SvE13 - C13; // S02
   SPK[2][0] = SPK[0][2];
}
