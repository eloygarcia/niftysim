// =========================================================================
// File:       tledContactUSProbe.h
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

#ifndef tledContactUSProbe_H
#define tledContactUSProbe_H

#ifdef _GPU_
#include "tledCUDAHelpers.h"
#endif

#include <vector>

#ifdef _GPU_
struct tledGPUContactUSProbe
{
   float4 Origin; // Contains: .x=Origin0, .y=Origin1, .z=Origin2, .w=Radius
   float4 Axis; // Contains: .x=Axis0, .y=Axis1, .z=Axis2, .w=Length
   int* SlaveMask;
};
#endif // _GPU_

/**
 * \ingroup contact 	 
 */
class tledContactUSProbe
{
public:
   tledContactUSProbe() {;}
   tledContactUSProbe(std::vector<float> orig, std::vector<float> axis, float R, float L, std::vector<int> slvs, std::vector<float> origdisp, float radchng, int NumNodes);
   ~tledContactUSProbe() {;}
   // Slaves
   std::vector<int>* GetSlvs() {return &slvNodes;}
   // Origin
   float* GetStartOrigin() {return Origin0;}
   std::vector<float> GetStartOriginV();
   float* GetCurrentOrigin() {return Origin;}
   void SetStartOrigin(std::vector<float> orig);
   // Axis
   float* GetStartAxis() {return Axis0;}
   std::vector<float> GetStartAxisV();
   float* GetCurrentAxis() {return Axis;}
   void SetStartAxis(std::vector<float> axis);
   // Radius
   float GetStartRadius() {return Radius0;}
   float GetCurrentRadius() {return Radius;}
   void SetStartRadius(float r) {Radius0 = r;}
   // Length
   float GetStartLength() {return Length0;}
   float GetCurrentLength() {return Length;}
   void SetStartLength(float l) {Length0 = l;}
   // OriginDisp
   float* GetOriginDisp() {return OriginDisp;}
   std::vector<float> GetOriginDispV();
   void SetOriginDisp(std::vector<float> disp);
   // RadiusChng
   float GetRadiusChng() {return RadiusChng;}
   void SetRadiusChng(float dr) {RadiusChng = dr;}
   
   // Update probe variables for relative simulation time TR
   void Update(double TR);
   
   // Return the allocated device variable pointer
#ifdef _GPU_
  /**
   * \name On-Device Memory
   * @{ 
   */
public:
  tledGPUContactUSProbe* GetGPUContactUSProbe(void) {return d_Prb;}

   /* Deallocation of GPU memory not automatically performed on destruction, must be done explicitly with this static member function. */
   static void ReleaseGPUMemory(tledGPUContactUSProbe *dp_prb);
   /** @} */
#endif // _GPU_
   
private:
   float Origin0[3];
   float Origin[3];
   float Axis0[3];
   float Axis[3];
   float Radius0;
   float Radius;
   float Length0;
   float Length;
   std::vector<int> slvNodes;
   float OriginDisp[3];
   float RadiusChng;
   
#ifdef _GPU_
   tledGPUContactUSProbe* d_Prb;
#endif // _GPU_
};

#endif // tledContactUSProbe_H

