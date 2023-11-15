// =========================================================================
// File:       tledContactPlate.h
// Purpose:    Create instance of a rigid rectangular plate for contact modelling
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

#ifndef tledContactPlate_H
#define tledContactPlate_H

#ifdef _GPU_
#include "tledCUDAHelpers.h"
#endif

#include <vector>

#ifdef _GPU_
struct tledGPUContactPlate
{
   float4 a;
   float4 b;
   float4 c;
   int* SlaveMask;
};
#endif // _GPU_

/**
 * \ingroup contact 	 
 */
class tledContactPlate
{
public:
   tledContactPlate() {;}
   tledContactPlate(std::vector<float> A, std::vector<float> B, std::vector<float> C, std::vector<int> slvs, std::vector<float> disp, int NumNodes);
   ~tledContactPlate() {;}
   // Slaves
   std::vector<int>* GetSlvs() {return &slvNodes;}
   // Corner A
   float* GetStartCrnrA() {return a0;}
   std::vector<float> GetStartCrnrAV();
   float* GetCurrentCrnrA() {return a;}
   void SetStartCrnrA(std::vector<float> A);
   // Corner B
   float* GetStartCrnrB() {return b0;}
   std::vector<float> GetStartCrnrBV();
   float* GetCurrentCrnrB() {return b;}
   void SetStartCrnrB(std::vector<float> B);
   // Corner C
   float* GetStartCrnrC() {return c0;}
   std::vector<float> GetStartCrnrCV();
   float* GetCurrentCrnrC() {return c;}
   void SetStartCrnrC(std::vector<float> C);
   // Disp
   float* GetDisp() {return Disp;}
   std::vector<float> GetDispV();
   void SetDisp(std::vector<float> disp);
   
   // Update plate variables for relative simulation time TR
   void Update(double TR);
   
#ifdef _GPU_
  /**
   * \name On-Device Memory
   * @{ 
   */
public:
   /** Return the allocated device variable pointer */
   tledGPUContactPlate* GetGPUContactPlate(void) {return d_Plt;}

   /* Deallocation of GPU memory not automatically performed on destruction, must be done explicitly with this static member function. */
   static void ReleaseGPUMemory(tledGPUContactPlate *dp_plt);
  /** @} */
#endif // _GPU_
   
private:
   // The rectangular plate is defined by 3 corner locations a,b,c in 3D space:
   //
   //               c---------d
   // z  y         /         /
   // | /         /         /
   // |/         /         /
   //  ---> x   a---------b
   //
   // The position of d is inferred. These corners also form a local coordinate system with origin at point a and
   // x- and y-axes in the direction of ab and ac, respectively. The local z-axis (= ab x ac) defines the plate
   // normal vector.
   
   float a0[3];
   float b0[3];
   float c0[3];
   float a[3];
   float b[3];
   float c[3];
   std::vector<int> slvNodes;
   float Disp[3];
   
#ifdef _GPU_
   tledGPUContactPlate* d_Plt;
#endif // _GPU_
};

#endif // tledContactPlate_H

