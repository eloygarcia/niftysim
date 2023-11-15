// =========================================================================
// File:       tledElementH8.cpp
// Purpose:    8-node reduced integration hexahedron element with hourglass control
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    April 2007
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#include "tledElementH8.h"
#include "tledMatrixFunctions.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

tledElementH8::~tledElementH8()
{
	delete Mat;
}

tledElementH8::tledElementH8(tledMesh* Mesh, int ElNum, const char* MatType, float* MatParams, float Density, float kappa)
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
   ENum = ElNum;
   EInd = Mesh->GetElNodeInds(ElNum);
   vector<float> NCds;
   for (int i = 0; i < 8; i++)
   {
      NCds = Mesh->GetNodeCds(EInd[i]);
      for (int j = 0; j < 3; j++)
      {
         x[i][j] = NCds[j];
      }
   }
   
   ComputeDhDx();
      
   HGKappa = kappa;
   ComputeHGParams();
}

void tledElementH8::ComputeDhDx()
{
   const float a = 1./8;
   float DhDr[8][3] = {{-a, -a, -a},	// Shape function natural derivatives
		       {a, -a, -a},
		       {a, a, -a},
		       {-a, a, -a},
		       {-a, -a, a},
		       {a, -a, a},
		       {a, a, a},
		       {-a, a, a}};
   // Jacobian
   float J[3][3];
   MatMult83T83(DhDr,x,J);
   // Jacobian determinant
   MatDet33(J,&detJ);
   // Compute element volume
   ComputeVol();
   // Compute shape function global derivatives
   float invJ[3][3];
   MatInv33(J,invJ);
   MatMult8333T(DhDr,invJ,DhDx);
}

void tledElementH8::SetGeometry(tledMesh* Mesh)
{
   vector<float> NCds;
   for (int i = 0; i < 8; i++)
   {
      NCds = Mesh->GetNodeCds(EInd[i]);
      for (int j = 0; j < 3; j++)
      {
         x[i][j] = NCds[j];
      }
   }

   ComputeDhDx();
   ComputeHGParams();
}

void tledElementH8::ComputeElementMass(float *M)
{
   int i;
   float Mass = Mat->GetDensity()*Vol/8;
   for (i = 0; i < 8; i++)
   {
      M[EInd[i]*3] += Mass;
      M[EInd[i]*3+1] += Mass;
      M[EInd[i]*3+2] += Mass;
   }
}

void tledElementH8::ComputeElementForces(const float* U, float* F)
{
   for (int i = 0; i < 8; i++)
   {
      int j = EInd[i]*3;
      u[i][0] = U[j+0]; u[i][1] = U[j+1]; u[i][2] = U[j+2];
   }
   // Displacement derivs
   MatMult83T83(u,DhDx,X);
   // Deformation gradient: X = DuDx + I
   X[0][0] += 1; X[1][1] += 1; X[2][2] += 1;	// X is now def. grad.
   // 2nd Piola-Kirchhoff stress
   Mat->ComputeSPKStress(X,SPK);
   // Compute element nodal forces: Fe = 8*detJ*X*S*dh'
   float temp[3][3];
   MatMult3333(X,SPK,temp);
   MatMultScalar(&temp[0][0],3,3,8*detJ,&temp[0][0]);
   MatMult3383T(temp,DhDx,Fe);

   // Compute HG forces
   ComputeHGForces();

   // Add element forces to global forces
   for (int i = 0; i < 8; i++)
   {
      int  j = EInd[i]*3;
      F[j] += Fe[0][i]; F[j+1] += Fe[1][i]; F[j+2] += Fe[2][i];
   }
}

void tledElementH8::ComputeHGParams()
{
   int i,j;
   float a = 0;
   for (i = 0; i < 8; i++)
   {
      for (j = 0; j < 3; j++)
      {
            a += DhDx[i][j]*DhDx[i][j];
      }
   }
   float L;
   float M;
   Mat->GetHGLameParams(&L,&M);
   float k = HGKappa*Vol*(L+2*M)*a/8;

   float Gamma[8][4] = {{1,1,1,-1},
                        {-1,1,-1,1},
                        {1,-1,-1,-1},
                        {-1,-1,1,1},
                        {1,-1,-1,1},
                        {-1,-1,1,-1},
                        {1,1,1,1},
                        {-1,1,-1,-1}};
   float A[8][8];
   MatMult8383T(DhDx,x,A);
   float gamma[8][4];
   MatMult8884(A,Gamma,gamma);
   MatSubtract(&Gamma[0][0],&gamma[0][0],8,4,&gamma[0][0]);

   MatMult8484T(gamma,gamma,HG);
   MatMultScalar(&HG[0][0],8,8,k,&HG[0][0]);
}

void tledElementH8::ComputeHGForces()
{
   MatMult8883(HG,u,FHG);
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 8; j++)
      {
         Fe[i][j] += FHG[j][i];
      }
   }
}

void tledElementH8::ComputeVol()
{
   // Calc CIJK first
   float C[8][8][8];
   ComputeCIJK(C);
   // Calc volume
   Vol = 0;
   int I,J,K;
   for (I = 0; I < 8; I++)
   {
      for (J = 0; J < 8; J++)
      {
         for (K = 0; K < 8; K++)
         {
            Vol += x[I][0]*x[J][1]*x[K][2]*C[I][J][K];
         }
      }
   }
}

void tledElementH8::ComputeDhDx_BmatMeth()
{
   // Calc B matrix
   float B[8][3];
   ComputeBmat(B);
   for (int i = 0; i < 8; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         DhDx[i][j] = B[i][j]/Vol;
      }
   }
}

void tledElementH8::ComputeBmat(float B[8][3])
{
   // Calc CIJK first
   float C[8][8][8];
   ComputeCIJK(C);
   // Calc B
   memset(B,0,sizeof(float)*8*3);
   for (int I = 0; I < 8; I++)
   {
      for (int J = 0; J < 8; J++)
      {
         for (int K = 0; K < 8; K++)
         {
            B[I][0] += x[J][1]*x[K][2]*C[I][J][K];
            B[I][1] += x[J][2]*x[K][0]*C[I][J][K];
            B[I][2] += x[J][0]*x[K][1]*C[I][J][K];
         }
      }
   }
}

void tledElementH8::ComputeCIJK(float C[8][8][8])
{
   float a = (float)(1./12);
   float Ctemp[8*8*8] = 
   {0,0,0,0,0,0,0,0,
   0,0,-a,-a,a,a,0,0,
   0,a,0,-a,0,0,0,0,
   0,a,a,0,-a,0,0,-a,
   0,-a,0,a,0,-a,0,a,
   0,-a,0,0,a,0,0,0,
   0,0,0,0,0,0,0,0,
   0,0,0,a,-a,0,0,0,

   0,0,a,a,-a,-a,0,0,
   0,0,0,0,0,0,0,0,
   -a,0,0,-a,0,a,a,0,
   -a,0,a,0,0,0,0,0,
   a,0,0,0,0,-a,0,0,
   a,0,-a,0,a,0,-a,0,
   0,0,-a,0,0,a,0,0,
   0,0,0,0,0,0,0,0,

   0,-a,0,a,0,0,0,0,
   a,0,0,a,0,-a,-a,0,
   0,0,0,0,0,0,0,0,
   -a,-a,0,0,0,0,a,a,
   0,0,0,0,0,0,0,0,
   0,a,0,0,0,0,-a,0,
   0,a,0,-a,0,a,0,-a,
   0,0,0,-a,0,0,a,0,

   0,-a,-a,0,a,0,0,a,
   a,0,-a,0,0,0,0,0,
   a,a,0,0,0,0,-a,-a,
   0,0,0,0,0,0,0,0,
   -a,0,0,0,0,0,0,a,
   0,0,0,0,0,0,0,0,
   0,0,a,0,0,0,0,-a,
   -a,0,a,0,-a,0,a,0,

   0,a,0,-a,0,a,0,-a,
   -a,0,0,0,0,a,0,0,
   0,0,0,0,0,0,0,0,
   a,0,0,0,0,0,0,-a,
   0,0,0,0,0,0,0,0,
   -a,-a,0,0,0,0,a,a,
   0,0,0,0,0,-a,0,a,
   a,0,0,a,0,-a,-a,0,

   0,a,0,0,-a,0,0,0,
   -a,0,a,0,-a,0,a,0,
   0,-a,0,0,0,0,a,0,
   0,0,0,0,0,0,0,0,
   a,a,0,0,0,0,-a,-a,
   0,0,0,0,0,0,0,0,
   0,-a,-a,0,a,0,0,a,
   0,0,0,0,a,0,-a,0,

   0,0,0,0,0,0,0,0,
   0,0,a,0,0,-a,0,0,
   0,-a,0,a,0,-a,0,a,
   0,0,-a,0,0,0,0,a,
   0,0,0,0,0,a,0,-a,
   0,a,a,0,-a,0,0,-a,
   0,0,0,0,0,0,0,0,
   0,0,-a,-a,a,a,0,0,

   0,0,0,-a,a,0,0,0,
   0,0,0,0,0,0,0,0,
   0,0,0,a,0,0,-a,0,
   a,0,-a,0,a,0,-a,0,
   -a,0,0,-a,0,a,a,0,
   0,0,0,0,-a,0,a,0,
   0,0,a,a,-a,-a,0,0,
   0,0,0,0,0,0,0,0};

   memcpy(C,Ctemp,sizeof(float)*8*8*8);
}

void tledElementH8::GetDhDx(float* Dh)
{
   for (int i = 0; i < 8; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         Dh[3*i + j] = DhDx[i][j];
      }
   }
}

void tledElementH8::GetStress(float* S)
{
   FillSv();
   memcpy(S,Sv,sizeof(float)*6);
}

void tledElementH8::GetStrain(float* E)
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

void tledElementH8::GetVMStress(float* SVM)
{
   FillSv();
   *SVM = sqrt( Sv[0]*Sv[0] + Sv[1]*Sv[1] + Sv[2]*Sv[2]
               - Sv[0]*Sv[1] - Sv[1]*Sv[2] - Sv[0]*Sv[2]
               + 3*(Sv[3]*Sv[3] + Sv[4]*Sv[4] + Sv[5]*Sv[5]) );
}

void tledElementH8::GetVMStrain(float* EVM)
{
   float E[6];
   this->GetStrain(E);
   float L,M,P;	// Lambda, Mu, PoissonRatio
   Mat->GetHGLameParams(&L,&M);
   P = L/(2*(L+M));
   *EVM = sqrt( E[0]*E[0] + E[1]*E[1] + E[2]*E[2]
               - E[0]*E[1] - E[1]*E[2] - E[0]*E[2]
               + 3*(E[3]*E[3] + E[4]*E[4] + E[5]*E[5]) )/(1+P);
}

void tledElementH8::FillSv()
{
   Sv[0] = SPK[0][0];
   Sv[1] = SPK[1][1];
   Sv[2] = SPK[2][2];
   Sv[3] = SPK[0][1];
   Sv[4] = SPK[1][2];
   Sv[5] = SPK[0][2];
   
//   cout << "SPK = " << Sv[0] << " " << Sv[1] << " " << Sv[2] << " " << Sv[3] << " " << Sv[4] << " " << Sv[5] << endl;
}

float tledElementH8::GetStrainEnergy(void)
{
   FillSv();
   float Ev[6];
   this->GetStrain(Ev);
   Ev[3] *= 2; Ev[4] *= 2; Ev[5] *= 2; // Shear components must be doubled
   // Strain energy: e = (1/2) * int_V {Dot(Ev,Sv)} dV = (1/2) * 8 * detJ * Dot(Ev,Sv) = 4 * detJ * Dot(Ev,Sv)
   return 4*detJ*Dot(Ev,6,Sv);
}

void tledElementH8::ComputeStress(float* U)
{
   for (int i = 0; i < 8; i++)
   {
      int j = EInd[i]*3;
      u[i][0] = U[j+0]; u[i][1] = U[j+1]; u[i][2] = U[j+2];
   }
   // Displacement derivs
   MatMult83T83(u,DhDx,X);
   // Deformation gradient: X = DuDx + I
   X[0][0] += 1; X[1][1] += 1; X[2][2] += 1;	// X is now def. grad.
   // 2nd Piola-Kirchhoff stress
   Mat->ComputeSPKStress(X,SPK);
}

void tledElementH8::SetMaterialParams(vector<float> params)
{
   Mat->SetParams(params);
   ComputeHGParams();
}
