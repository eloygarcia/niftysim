// =========================================================================
// File:       tledContactUSProbe.cu
// Purpose:    Create instance of a rigid ultrasound probe for contact
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    August 2010
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifdef _GPU_

#include "tledContactUSProbe.h"

#include <iostream>

using namespace std;

tledContactUSProbe::tledContactUSProbe(vector<float> orig, vector<float> axis, float R, float L, vector<int> slvs, vector<float> origdisp, float radchng, int NumNodes)
{
   for (int i = 0; i < 3; i++)
   {
      Origin0[i] = orig[i];
      Axis0[i] = axis[i];
      OriginDisp[i] = origdisp[i];
   }
   memcpy(Origin,Origin0,sizeof(float)*3);
   memcpy(Axis,Axis0,sizeof(float)*3);
   Radius0 = Radius = R;
   Length0 = Length = L;
   slvNodes = slvs;
   RadiusChng = radchng;
   
   // Allocate device variable
   
   tledGPUContactUSProbe l_Prb;
   l_Prb.Origin.x = Origin[0];
   l_Prb.Origin.y = Origin[1];
   l_Prb.Origin.z = Origin[2];
   l_Prb.Origin.w = Radius;
   l_Prb.Axis.x = Axis[0];
   l_Prb.Axis.y = Origin[1];
   l_Prb.Axis.z = Origin[2];
   l_Prb.Axis.w = Radius;
   cudaMalloc((void**)&(l_Prb.SlaveMask),sizeof(int)*NumNodes);
   int* l_SlaveMask = new int[NumNodes];
   memset(l_SlaveMask,0,sizeof(int)*NumNodes);
   for (int i = 0; i < (int)slvNodes.size(); i++)
      l_SlaveMask[slvNodes[i]] = 1;
   cudaMemcpy(l_Prb.SlaveMask,l_SlaveMask,sizeof(int)*NumNodes,cudaMemcpyHostToDevice);
   
   cudaMalloc((void**)&d_Prb,sizeof(tledGPUContactUSProbe));
   
   cudaMemcpy(d_Prb,&l_Prb,sizeof(tledGPUContactUSProbe),cudaMemcpyHostToDevice);
   delete l_SlaveMask;
}

void tledContactUSProbe::ReleaseGPUMemory(tledGPUContactUSProbe *dp_prb) {
  int *dp_slaveMask = NULL;

  tledCUDAHelpers::CopyFromDevice<int*>(&dp_slaveMask, &dp_prb->SlaveMask);
  tledCheckCUDAErrors(cudaFree(dp_slaveMask));
  tledCheckCUDAErrors(cudaFree(dp_prb));
}

void tledContactUSProbe::Update(double TR)
{
   double TR2 = TR*TR;
   float Amp = (float)( TR2*(10*TR - TR2*(15 - 6*TR)) );
   // Update origin
   for (int j = 0; j < 3; j++)
      Origin[j] = Origin0[j] + OriginDisp[j]*Amp;
   // Update radius
   Radius = Radius0 + RadiusChng*Amp;
   // Update device variables
   float4 l_Origin;
   l_Origin.x = Origin[0];
   l_Origin.y = Origin[1];
   l_Origin.z = Origin[2];
   l_Origin.w = Radius;
   cudaMemcpy(&(d_Prb->Origin),&l_Origin,sizeof(float4),cudaMemcpyHostToDevice);
   float4 l_Axis;
   l_Axis.x = Axis[0];
   l_Axis.y = Axis[1];
   l_Axis.z = Axis[2];
   l_Axis.w = Length;
   cudaMemcpy(&(d_Prb->Axis),&l_Axis,sizeof(float4),cudaMemcpyHostToDevice);
}

vector<float> tledContactUSProbe::GetStartOriginV()
{
   vector<float> orig(3);
   orig[0] = Origin0[0]; orig[1] = Origin0[1]; orig[2] = Origin0[2];
   return orig;
}

vector<float> tledContactUSProbe::GetStartAxisV()
{
   vector<float> axis(3);
   axis[0] = Axis0[0]; axis[1] = Axis0[1]; axis[2] = Axis0[2];
   return axis;
}

vector<float> tledContactUSProbe::GetOriginDispV()
{
   vector<float> disp(3);
   disp[0] = OriginDisp[0]; disp[1] = OriginDisp[1]; disp[2] = OriginDisp[2];
   return disp;
}

void tledContactUSProbe::SetOriginDisp(vector<float> disp)
{
   if (disp.size() != 3)
   {
      cerr << "!!! Warning: invalid vector passed" << endl;
      return;
   }
   for (int i = 0; i < 3; i++)
      OriginDisp[i] = disp[i];
}

void tledContactUSProbe::SetStartOrigin(vector<float> orig)
{
   if (orig.size() != 3)
   {
      cerr << "!!! Warning: invalid vector passed" << endl;
      return;
   }
   for (int i = 0; i < 3; i++)
      Origin0[i] = orig[i];
}

void tledContactUSProbe::SetStartAxis(vector<float> axis)
{
   if (axis.size() != 3)
   {
      cerr << "!!! Warning: invalid vector passed" << endl;
      return;
   }
   for (int i = 0; i < 3; i++)
      Axis0[i] = axis[i];
}

#endif // _GPU_

