// =========================================================================
// File:       tledContactCylinder.h
// Purpose:    Create instance of a rigid cylinder for contact
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    July 2010
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifndef tledContactCylinder_H
#define tledContactCylinder_H

#ifdef _GPU_
#include "tledCUDAHelpers.h"
#endif

#include <vector>
#include <cstring>

#ifdef _GPU_

struct tledGPUContactCylinder
{
   float4 Origin; // Contains: .x=Origin0, .y=Origin1, .z=Origin2, .w=Radius
   float4 Axis; // Contains: .x=Axis0, .y=Axis1, .z=Axis2, .w=Length
   int* SlaveMask;
};
#endif // _GPU_

/**
 * \ingroup contact 	 
 */
class tledContactCylinder
{
public:
   tledContactCylinder() {;}
   tledContactCylinder(std::vector<float> orig, std::vector<float> axis, float R, float L, std::vector<int> slvs, std::vector<float> origdisp, float radchng, int NumNodes);
   ~tledContactCylinder() {}
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
   // Length
   float GetStartLength() {return Length0;}
   float GetCurrentLength() {return Length;}
   // OriginDisp
   float* GetOriginDisp() {return OriginDisp;}
   std::vector<float> GetOriginDispV();
   void SetOriginDisp(std::vector<float> disp);
   // RadiusChng
   float GetRadiusChng() {return RadiusChng;}
   void SetRadiusChng(float dr) {RadiusChng = dr;}
   void SetStartRadius(float r) {Radius0 = r;}
   // Length
   void SetStartLength(float l) {Length0 = l;}
   
   // Update cylinder variables for relative simulation time TR
   void Update(double TR);
   
#ifdef _GPU_
  /**
   * \name On-Device Memory
   * @{ 
   */
public:
   /* Return the allocated device variable pointer */
   tledGPUContactCylinder* GetGPUContactCylinder(void) {return d_Cyl;}

   /* Deallocation of GPU memory not automatically performed on destruction, must be done explicitly with this static member function. */
   static void ReleaseGPUMemory(tledGPUContactCylinder *dp_cyl);
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
   tledGPUContactCylinder* d_Cyl;
#endif // _GPU_
};


#endif // tledContactCylinder_H

