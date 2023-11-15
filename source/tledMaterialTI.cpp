// =========================================================================
// File:       tledMaterialTI.cpp
// Purpose:    Material class - Transversely isotropic
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

#include "tledMaterialTI.h"
#include "tledMatrixFunctions.h"

#include <iostream>
#include <cmath>

using namespace std;

tledMaterialTI::tledMaterialTI()
{
}

tledMaterialTI::tledMaterialTI(float* MatParams, float Density)
{
   this->SetDensity(Density);

   // MatParams[0] = M, MatParams[1] = K, MatParams[2] = E, MatParams[3->5] = a
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
}

tledMaterialTI::~tledMaterialTI()
{
}

void tledMaterialTI::GetHGLameParams(float* Lambda, float* Mu)
{
   *Lambda = K - 2*M/3;
   *Mu = M;
}

void tledMaterialTI::SetParams(vector<float> params)
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
}

void tledMaterialTI::ComputeSPKStress(float X[3][3], float *SPK)
{
   float J,invdetC,x1,x2,x3,x4,x5;
   float XT11,XT12,XT13,XT21,XT22,XT23,XT31,XT32,XT33;
   float C11,C12,C13,C22,C23,C33;
   float Ci11,Ci12,Ci13,Ci22,Ci23,Ci33;
   float Si11,Si12,Si13,Si22,Si23,Si33;
   float Sv11,Sv12,Sv13,Sv22,Sv23,Sv33;

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
   x5 = K*J*(J-1);
   J = (float)pow((double)J,-(double)2/(double)3);	// J23 = J^(-2/3)
   x1 = J*M;
   x2 = J*(A11*C11+A22*C22+A33*C33+2*A12*C12+2*A23*C23+2*A13*C13) - 1; // Bracketed term is I4 = A:C
   x3 = J*E*x2;
   x4 = -(E*x2*(x2+1) + x1*(C11+C22+C33))/3;

   // C inverse
   invdetC = 1/(C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23));
   Ci11 = (C22*C33 - C23*C23)*invdetC;
   Ci12 = (C13*C23 - C12*C33)*invdetC;
   Ci13 = (C12*C23 - C13*C22)*invdetC;
   Ci22 = (C11*C33 - C13*C13)*invdetC;
   Ci23 = (C12*C13 - C11*C23)*invdetC;
   Ci33 = (C11*C22 - C12*C12)*invdetC;

   // Isochoric stresses
   Si11 = x3*A11 + x4*Ci11 + x1;
   Si22 = x3*A22 + x4*Ci22 + x1;
   Si33 = x3*A33 + x4*Ci33 + x1;
   Si12 = x3*A12 + x4*Ci12;
   Si23 = x3*A23 + x4*Ci23;
   Si13 = x3*A13 + x4*Ci13;

   // Volumetric stresses
   Sv11 = x5*Ci11;
   Sv22 = x5*Ci22;
   Sv33 = x5*Ci33;
   Sv12 = x5*Ci12;
   Sv23 = x5*Ci23;
   Sv13 = x5*Ci13;

   // Total stresses
   SPK[0] = Si11 + Sv11;
   SPK[1] = Si22 + Sv22;
   SPK[2] = Si33 + Sv33;
   SPK[3] = Si12 + Sv12;
   SPK[4] = Si23 + Sv23;
   SPK[5] = Si13 + Sv13;
}

void tledMaterialTI::ComputeSPKStress(float X[3][3], float SPK[3][3])
{
   float J,invdetC,x1,x2,x3,x4,x5;
   float XT11,XT12,XT13,XT21,XT22,XT23,XT31,XT32,XT33;
   float C11,C12,C13,C22,C23,C33;
   float Ci11,Ci12,Ci13,Ci22,Ci23,Ci33;
   float Si11,Si12,Si13,Si22,Si23,Si33;
   float Sv11,Sv12,Sv13,Sv22,Sv23,Sv33;

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
   x5 = K*J*(J-1);
   J = (float)pow((double)J,-(double)2/(double)3);	// J23 = J^(-2/3)
   x1 = J*M;
   x2 = J*(A11*C11+A22*C22+A33*C33+2*A12*C12+2*A23*C23+2*A13*C13) - 1; // Bracketed term is I4 = A:C
   x3 = J*E*x2;
   x4 = -(E*x2*(x2+1) + x1*(C11+C22+C33))/3;

   // C inverse
   invdetC = 1/(C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23));
   Ci11 = (C22*C33 - C23*C23)*invdetC;
   Ci12 = (C13*C23 - C12*C33)*invdetC;
   Ci13 = (C12*C23 - C13*C22)*invdetC;
   Ci22 = (C11*C33 - C13*C13)*invdetC;
   Ci23 = (C12*C13 - C11*C23)*invdetC;
   Ci33 = (C11*C22 - C12*C12)*invdetC;

   // Isochoric stresses
   Si11 = x3*A11 + x4*Ci11 + x1;
   Si22 = x3*A22 + x4*Ci22 + x1;
   Si33 = x3*A33 + x4*Ci33 + x1;
   Si12 = x3*A12 + x4*Ci12;
   Si23 = x3*A23 + x4*Ci23;
   Si13 = x3*A13 + x4*Ci13;

   // Volumetric stresses
   Sv11 = x5*Ci11;
   Sv22 = x5*Ci22;
   Sv33 = x5*Ci33;
   Sv12 = x5*Ci12;
   Sv23 = x5*Ci23;
   Sv13 = x5*Ci13;

   // Total stress
   SPK[0][0] = Si11 + Sv11; // S00
   SPK[1][1] = Si22 + Sv22; // S11
   SPK[2][2] = Si33 + Sv33; // S22
   SPK[0][1] = Si12 + Sv12; // S01
   SPK[1][0] = SPK[0][1];
   SPK[1][2] = Si23 + Sv23; // S12
   SPK[2][1] = SPK[1][2];
   SPK[0][2] = Si13 + Sv13; // S02
   SPK[2][0] = SPK[0][2];
}
