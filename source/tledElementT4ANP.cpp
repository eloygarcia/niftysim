// =========================================================================
// File:       tledElementT4ANP.cpp
// Purpose:    4-node tetrahedron element which uses improved average nodal pressure formulation
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    March 2009
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#include "tledElementT4ANP.h"
#include "tledMatrixFunctions.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

tledElementT4ANP::~tledElementT4ANP()
{
   delete Mat;
}

tledElementT4ANP::tledElementT4ANP(tledMesh* Mesh, int ElNum, const char* MatType, float* MatParams, float Density)
{
   // Initialise element characteristics
   if (!strcmp(MatType,"LE"))
   {
      Mat = new tledMaterialLE(MatParams,Density);
   }
   else if (!strcmp(MatType,"NH"))
   {
      Mat = new tledMaterialNH(MatParams,Density);
   }
   else if (!strcmp(MatType,"AB"))
   {
      Mat = new tledMaterialAB(MatParams,Density);
   }
   else if (!strcmp(MatType,"PY"))
   {
      Mat = new tledMaterialPY(MatParams,Density);
   }
   else if (!strcmp(MatType,"NHV"))
   {
      Mat = new tledMaterialNHV(MatParams,Density);
   }
   else if (!strcmp(MatType,"TI"))
   {
      Mat = new tledMaterialTI(MatParams,Density);
   }
   else if (!strcmp(MatType,"TIV"))
   {
      Mat = new tledMaterialTIV(MatParams,Density);
   }
   NDOF = 3;
   EInd = Mesh->GetElNodeInds(ElNum);
   vector<float> NCds;
   for (int i = 0; i < 4; i++)
   {
      NCds = Mesh->GetNodeCds(EInd[i]);
      for (int j = 0; j < 3; j++)
         x[i][j] = NCds[j];
   }

   ComputeDhDx();
}

void tledElementT4ANP::ComputeDhDx()
{
  float DhDr[4][3] = {{-1, -1, -1},	// Shape function natural derivatives
		      {1, 0, 0},
		      {0, 1, 0},
		      {0, 0, 1}};
   float Jac[3][3];	// Jacobian
   MatMult43T43(DhDr,x,Jac);
   // Compute element volume - note alternative, but equivalent formula
   float detJ;
   MatDet33(Jac,&detJ);
   Vol = fabs(detJ/6);
   // Compute shape function global derivatives
   float invJ[3][3];
   MatInv33(Jac,invJ);
   MatMult4333T(DhDr,invJ,DhDx);
}

void tledElementT4ANP::SetGeometry(tledMesh* Mesh)
{
   vector<float> NCds;
   for (int i = 0; i < 4; i++)
   {
      NCds = Mesh->GetNodeCds(EInd[i]);
      for (int j = 0; j < 3; j++)
      {
         x[i][j] = NCds[j];
      }
   }
   
   ComputeDhDx();
}

void tledElementT4ANP::ComputeElementMass(float *M)
{
   int i;
   float Mass = Mat->GetDensity()*Vol/4;
   for (i = 0; i < 4; i++)
   {
      M[EInd[i]*3] += Mass;
      M[EInd[i]*3+1] += Mass;
      M[EInd[i]*3+2] += Mass;
   }
}

void tledElementT4ANP::ComputeElementPressure(const float* U, float* Pa)
{
   for (int i = 0; i < 4; i++)
   {
      int j = EInd[i]*3;
      u[i][0] = U[j+0]; u[i][1] = U[j+1]; u[i][2] = U[j+2];
   }
   // Displacement derivs
   MatMult43T43(u,DhDx,X);
   // Deformation gradient: X = DuDx + I
   X[0][0] += 1; X[1][1] += 1; X[2][2] += 1;	// X is now def. grad.
   // Element Jacobian
   MatDet33(X,&J);
   // Element pressure
   float kappa;
   Mat->GetKappa(&kappa);
   float P = kappa*(J-1);
   // Add to nodal sums
   for (int i = 0; i < 4; i++)
      Pa[EInd[i]] += P*Vol;
}

void tledElementT4ANP::ComputeModifiedElementForces(const float* Pa, const float* Va, float* F)
{
   // Modified pressure
   float Pm = 0;
   for (int i = 0; i < 4; i++)
      Pm += Pa[EInd[i]]/Va[EInd[i]];
   Pm /= 16;
   // Modified Jacobian
   float kappa;
   Mat->GetKappa(&kappa);
   float Jm = Pm/kappa + 1;
   // Modified def'n gradient
   Jm = (float)std::pow((double)Jm,(double)1/(double)3);	// Jm = Jm^1/3
   J = (float)std::pow((double)J,-(double)1/(double)3);	// J = J^-1/3
   MatMultScalar(&X[0][0],3,3,Jm*J,&X[0][0]);
   // 2nd Piola-Kirchhoff stress
   Mat->ComputeSPKStress(X,SPK);
   // Compute element nodal forces: Fe = V*X*S*dh'
   float temp[3][3];
   MatMult3333(X,SPK,temp);
   MatMultScalar(&temp[0][0],3,3,Vol,&temp[0][0]);
   MatMult3343T(temp,DhDx,Fe);
   // Add element forces to global forces
   for (int i = 0; i < 4; i++)
   {
      int  j = EInd[i]*3;
      F[j] += Fe[0][i]; F[j+1] += Fe[1][i]; F[j+2] += Fe[2][i];
   }
}

void tledElementT4ANP::GetDhDx(float* Dh)
{
   for (int i = 0; i < 4; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         Dh[3*i + j] = DhDx[i][j];
      }
   }
}

void tledElementT4ANP::GetStress(float* S)
{
   FillSv();
   memcpy(S,Sv,sizeof(float)*6);
}

void tledElementT4ANP::GetStrain(float* E)
{
   float XT11,XT12,XT13,XT21,XT22,XT23,XT31,XT32,XT33;
   float C11,C12,C13,C22,C23,C33;

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

   // Green-Lagrange strain
   E[0] = (C11 - 1)/2;
   E[1] = (C22 - 1)/2;
   E[2] = (C33 - 1)/2;
   E[3] = C12/2;
   E[4] = C23/2;
   E[5] = C13/2;
}

void tledElementT4ANP::GetVMStress(float* SVM)
{
   FillSv();
   *SVM = std::sqrt( Sv[0]*Sv[0] + Sv[1]*Sv[1] + Sv[2]*Sv[2]
		     - Sv[0]*Sv[1] - Sv[1]*Sv[2] - Sv[0]*Sv[2]
		     + 3*(Sv[3]*Sv[3] + Sv[4]*Sv[4] + Sv[5]*Sv[5]) );
}

void tledElementT4ANP::GetVMStrain(float* EVM)
{
   float E[6];
   GetStrain(E);
   float L,M,P;	// Lambda, Mu, PoissonRatio
   Mat->GetHGLameParams(&L,&M);
   P = L/(2*(L+M));
   *EVM = std::sqrt( E[0]*E[0] + E[1]*E[1] + E[2]*E[2]
		     - E[0]*E[1] - E[1]*E[2] - E[0]*E[2]
		     + 3*(E[3]*E[3] + E[4]*E[4] + E[5]*E[5]) )/(1+P);
}

void tledElementT4ANP::FillSv()
{
   Sv[0] = SPK[0][0];
   Sv[1] = SPK[1][1];
   Sv[2] = SPK[2][2];
   Sv[3] = SPK[0][1];
   Sv[4] = SPK[1][2];
   Sv[5] = SPK[0][2];
}

float tledElementT4ANP::GetStrainEnergy(void)
{
   FillSv();
   float Ev[6];
   this->GetStrain(Ev);
   Ev[3] *= 2; Ev[4] *= 2; Ev[5] *= 2; // Shear components must be doubled
   // Strain energy: e = (1/2) * int_V {Dot(Ev,Sv)} dV = (1/2) * V * Dot(Ev,Sv)
   return 0.5f*Vol*Dot(Ev,6,Sv);
}

void tledElementT4ANP::ComputeModifiedStress(float* Pa, float* Va)
{
   // Modified pressure
   float Pm = 0;
   for (int i = 0; i < 4; i++)
      Pm += Pa[EInd[i]]/Va[EInd[i]];
   Pm /= 16;
   // Modified Jacobian
   float kappa;
   Mat->GetKappa(&kappa);
   float Jm = Pm/kappa + 1;
   // Modified def'n gradient
   Jm = (float)std::pow((double)Jm,(double)1/(double)3);	// Jm = Jm^1/3
   J = (float)std::pow((double)J,-(double)1/(double)3);	// J = J^-1/3
   MatMultScalar(&X[0][0],3,3,Jm*J,&X[0][0]);
   // 2nd Piola-Kirchhoff stress
   Mat->ComputeSPKStress(X,SPK);
}

void tledElementT4ANP::SetMaterialParams(vector<float> params)
{
   Mat->SetParams(params);
}
