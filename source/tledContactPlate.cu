// =========================================================================
// File:       tledContactPlate.cu
// Purpose:    Create instance of a rigid rectangular sheet for contact modelling
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    November 2010
//
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
//
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifdef _GPU_

#include "tledContactPlate.h"

#include <iostream>

using namespace std;

tledContactPlate::tledContactPlate(vector<float> A, vector<float> B, vector<float> C, vector<int> slvs, vector<float> disp, int NumNodes)
{
   for (int i = 0; i < 3; i++)
   {
      a0[i] = A[i];
      a[i] = A[i];
      b0[i] = B[i];
      b[i] = B[i];
      c0[i] = C[i];
      c[i] = C[i];
      Disp[i] = disp[i];
   }
   slvNodes = slvs;
   
   // Allocate device variable
   
   tledGPUContactPlate l_Plt;
   l_Plt.a.x = a[0]; l_Plt.a.y = a[1]; l_Plt.a.y = a[2]; l_Plt.a.w = 0;
   l_Plt.b.x = b[0]; l_Plt.b.y = b[1]; l_Plt.b.y = b[2]; l_Plt.b.w = 0;
   l_Plt.c.x = c[0]; l_Plt.c.y = c[1]; l_Plt.c.y = c[2]; l_Plt.c.w = 0;
   cudaMalloc((void**)&(l_Plt.SlaveMask),sizeof(int)*NumNodes);
   int* l_SlaveMask = new int[NumNodes];
   memset(l_SlaveMask,0,sizeof(int)*NumNodes);
   for (int i = 0; i < (int)slvNodes.size(); i++)
      l_SlaveMask[slvNodes[i]] = 1;
   cudaMemcpy(l_Plt.SlaveMask,l_SlaveMask,sizeof(int)*NumNodes,cudaMemcpyHostToDevice);
   
   cudaMalloc((void**)&d_Plt,sizeof(tledGPUContactPlate));
   
   cudaMemcpy(d_Plt,&l_Plt,sizeof(tledGPUContactPlate),cudaMemcpyHostToDevice);
   delete l_SlaveMask;
}

void tledContactPlate::ReleaseGPUMemory(tledGPUContactPlate *dp_plt) {
  int *dp_slaveMask = NULL;

  tledCUDAHelpers::CopyFromDevice<int*>(&dp_slaveMask, &dp_plt->SlaveMask);
  tledCheckCUDAErrors(cudaFree(dp_slaveMask));
  tledCheckCUDAErrors(cudaFree(dp_plt));
}

void tledContactPlate::Update(double TR)
{
   double TR2 = TR*TR;
   float Amp = (float)( TR2*(10*TR - TR2*(15 - 6*TR)) );
   // Update corners
   for (int j = 0; j < 3; j++)
   {
      a[j] = a0[j] + Disp[j]*Amp;
      b[j] = b0[j] + Disp[j]*Amp;
      c[j] = c0[j] + Disp[j]*Amp;
   }
   // Update device variables
   float4 l_a;
   l_a.x = a[0]; l_a.y = a[1]; l_a.z = a[2];
   cudaMemcpy(&(d_Plt->a),&l_a,sizeof(float4),cudaMemcpyHostToDevice);
   float4 l_b;
   l_b.x = b[0]; l_b.y = b[1]; l_b.z = b[2];
   cudaMemcpy(&(d_Plt->b),&l_b,sizeof(float4),cudaMemcpyHostToDevice);
   float4 l_c;
   l_c.x = c[0]; l_c.y = c[1]; l_c.z = c[2];
   cudaMemcpy(&(d_Plt->c),&l_c,sizeof(float4),cudaMemcpyHostToDevice);
}

vector<float> tledContactPlate::GetStartCrnrAV()
{
   vector<float> A(3);
   A[0] = a0[0]; A[1] = a0[1]; A[2] = a0[2];
   return A;
}

vector<float> tledContactPlate::GetStartCrnrBV()
{
   vector<float> B(3);
   B[0] = b0[0]; B[1] = b0[1]; B[2] = b0[2];
   return B;
}

vector<float> tledContactPlate::GetStartCrnrCV()
{
   vector<float> C(3);
   C[0] = c0[0]; C[1] = c0[1]; C[2] = c0[2];
   return C;
}

vector<float> tledContactPlate::GetDispV()
{
   vector<float> disp(3);
   disp[0] = Disp[0]; disp[1] = Disp[1]; disp[2] = Disp[2];
   return disp;
}

void tledContactPlate::SetDisp(vector<float> disp)
{
   if (disp.size() != 3)
   {
      cerr << "!!! Warning: invalid vector passed" << endl;
      return;
   }
   for (int i = 0; i < 3; i++)
      Disp[i] = disp[i];
}

void tledContactPlate::SetStartCrnrA(vector<float> A)
{
   if (A.size() != 3)
   {
      cerr << "!!! Warning: invalid vector passed" << endl;
      return;
   }
   for (int i = 0; i < 3; i++)
      a0[i] = A[i];
}

void tledContactPlate::SetStartCrnrB(vector<float> B)
{
   if (B.size() != 3)
   {
      cerr << "!!! Warning: invalid vector passed" << endl;
      return;
   }
   for (int i = 0; i < 3; i++)
      b0[i] = B[i];
}

void tledContactPlate::SetStartCrnrC(vector<float> C)
{
   if (C.size() != 3)
   {
      cerr << "!!! Warning: invalid vector passed" << endl;
      return;
   }
   for (int i = 0; i < 3; i++)
      c0[i] = C[i];
}

#endif // _GPU_

