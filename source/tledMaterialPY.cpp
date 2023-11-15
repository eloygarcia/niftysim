// =========================================================================
// File:       tledMaterialPY.cpp
// Purpose:    Material class - Polynomial (N=2)
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    October 2010
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#include "tledMaterialPY.h"
#include "tledMatrixFunctions.h"

#include <iostream>
#include <cmath>

using namespace std;

tledMaterialPY::tledMaterialPY()
{
}

tledMaterialPY::tledMaterialPY(float* MatParams, float Density) : tledMaterial(Density)
{      
   c10 = MatParams[0];
   c01 = MatParams[1];
   c20 = MatParams[2];
   c02 = MatParams[3];
   c11 = MatParams[4];
   K = MatParams[5]; // Initial bulk modulus
   M = 2*(c10+c01); // Initial shear modulus
}

tledMaterialPY::~tledMaterialPY()
{
}

void tledMaterialPY::GetHGLameParams(float* Lambda, float* Mu)
{
   *Lambda = K - 2*M/3;
   *Mu = M;
}

void tledMaterialPY::SetParams(vector<float> params)
{
   c10 = params[0];
   c01 = params[1];
   c20 = params[2];
   c02 = params[3];
   c11 = params[4];
   K = params[5]; // Initial bulk modulus
   M = 2*(c10+c01); // Initial shear modulus
}

void tledMaterialPY::ComputeSPKStress(float X[3][3], float *SPK)
{
   //
   // EXPLANATION:
   //
   // S = Siso + Svol
   // where
   // Siso = 2*J^(-2/3)*(c10*A1 + c01*A2 + 2*C20*a1*A1 + 2*C02*a2*A2 + C11*(a2*A1+a1*A2))
   // with
   // a1 = I1b-3
   // a2 = I2b-3
   // A1 = I - I1*Ci/3
   // A2 = I1b*I - Cb - 2*J^(-2/3)*I2*Ci/3
   // and
   // Svol = K*J*(J-1)*Ci
   //
   
   // Transpose of deformation gradient
   float XT11 = X[0][0]; float XT12 = X[1][0]; float XT13 = X[2][0];
   float XT21 = X[0][1]; float XT22 = X[1][1]; float XT23 = X[2][1];
   float XT31 = X[0][2]; float XT32 = X[1][2]; float XT33 = X[2][2];
   
   // Right Cauchy-Green deformation tensor
   float C11 = XT11*XT11 + XT12*XT12 + XT13*XT13;
   float C12 = XT11*XT21 + XT12*XT22 + XT13*XT23;
   float C13 = XT11*XT31 + XT12*XT32 + XT13*XT33;
   float C22 = XT21*XT21 + XT22*XT22 + XT23*XT23;
   float C23 = XT21*XT31 + XT22*XT32 + XT23*XT33;
   float C33 = XT31*XT31 + XT32*XT32 + XT33*XT33;
   
   // Determinant of deformation gradient
   float J = XT11*(XT22*XT33 - XT32*XT23) + XT12*(XT23*XT31 - XT21*XT33) + XT13*(XT21*XT32 - XT22*XT31);
   float J23 = (float)pow((double)J,-(double)2/(double)3);
   // Determinant of C
   float detC = C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23);
   // Invariants
   float I1 = C11+C22+C33;
   float I2 = (float)(0.5*(I1*I1 - (C11*C11 + 2*C12*C12 + 2*C13*C13 + C22*C22 + 2*C23*C23 + C33*C33)));
   float I1b = J23*I1;
   float I2b = J23*J23*I2;
   
   // Some convenience variables
   float a1 = I1b-3;
   float a2 = I2b-3;
   float a3 = 2*J23*I2/3;
   float T1 = 2*J23*(c10 + c01*I1b + 2*c20*a1 + 2*c02*a2*I1b + c11*(a2 + a1*I1b));
   float T2 = -2*J23*J23*(c01 + 2*c02*a2 + c11*a1);
   float T3 = ( K*J*(J-1)
               -2*J23*(c10*I1/3 + c01*a3 + 2*c20*a1*I1/3 + 2*c02*a2*a3 + c11*(a2*I1/3 + a1*a3))
               )/detC;
   
   SPK[0] = T2*C11 + T3*(C22*C33-C23*C23) + T1;
   SPK[1] = T2*C22 + T3*(C11*C33-C13*C13) + T1;
   SPK[2] = T2*C33 + T3*(C11*C22-C12*C12) + T1;
   SPK[3] = T2*C12 + T3*(C13*C23-C12*C33);
   SPK[4] = T2*C23 + T3*(C12*C13-C11*C23);
   SPK[5] = T2*C13 + T3*(C12*C23-C13*C22);
}

void tledMaterialPY::ComputeSPKStress(float X[3][3], float SPK[3][3])
{
   // Transpose of deformation gradient
   float XT11 = X[0][0]; float XT12 = X[1][0]; float XT13 = X[2][0];
   float XT21 = X[0][1]; float XT22 = X[1][1]; float XT23 = X[2][1];
   float XT31 = X[0][2]; float XT32 = X[1][2]; float XT33 = X[2][2];
   
   // Right Cauchy-Green deformation tensor
   float C11 = XT11*XT11 + XT12*XT12 + XT13*XT13;
   float C12 = XT11*XT21 + XT12*XT22 + XT13*XT23;
   float C13 = XT11*XT31 + XT12*XT32 + XT13*XT33;
   float C22 = XT21*XT21 + XT22*XT22 + XT23*XT23;
   float C23 = XT21*XT31 + XT22*XT32 + XT23*XT33;
   float C33 = XT31*XT31 + XT32*XT32 + XT33*XT33;
   
   // Determinant of deformation gradient
   float J = XT11*(XT22*XT33 - XT32*XT23) + XT12*(XT23*XT31 - XT21*XT33) + XT13*(XT21*XT32 - XT22*XT31);
   float J23 = (float)pow((double)J,-(double)2/(double)3);
   // Determinant of C
   float detC = C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23);
   // Invariants
   float I1 = C11+C22+C33;
   float I2 = (float)(0.5*(I1*I1 - (C11*C11 + 2*C12*C12 + 2*C13*C13 + C22*C22 + 2*C23*C23 + C33*C33)));
   float I1b = J23*I1;
   float I2b = J23*J23*I2;

   // Some convenience variables
   float a1 = I1b-3;
   float a2 = I2b-3;
   float a3 = 2*J23*I2/3;
   float T1 = 2*J23*(c10 + c01*I1b + 2*c20*a1 + 2*c02*a2*I1b + c11*(a2 + a1*I1b));
   float T2 = -2*J23*J23*(c01 + 2*c02*a2 + c11*a1);
   float T3 = ( K*J*(J-1)
               -2*J23*(c10*I1/3 + c01*a3 + 2*c20*a1*I1/3 + 2*c02*a2*a3 + c11*(a2*I1/3 + a1*a3))
               )/detC;
   
   SPK[0][0] = T2*C11 + T3*(C22*C33-C23*C23) + T1;
   SPK[1][1] = T2*C22 + T3*(C11*C33-C13*C13) + T1;
   SPK[2][2] = T2*C33 + T3*(C11*C22-C12*C12) + T1;
   SPK[0][1] = T2*C12 + T3*(C13*C23-C12*C33);
   SPK[1][0] = SPK[0][1];
   SPK[1][2] = T2*C23 + T3*(C12*C13-C11*C23);
   SPK[2][1] = SPK[1][2];
   SPK[0][2] = T2*C13 + T3*(C12*C23-C13*C22);
   SPK[2][0] = SPK[0][2];
   
//   cout << "SPK = " << SPK[0][0] << " " << SPK[1][1] << " " << SPK[2][2] << endl;
}
