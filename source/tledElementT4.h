// =========================================================================
// File:       tledElementT4.h
// Purpose:    4-node tetrahedron element
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    March 2007
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledElementT4_H
#define tledElementT4_H

#include "tledElement.h"
#include "tledMesh.h"

class tledElementT4 : public tledElement
{
public:
   tledElementT4() {;}
   tledElementT4(tledMesh* Mesh, int ElNum, const char* MatType, float* MatParams, float Density);
   virtual ~tledElementT4();
   
   virtual void ComputeElementMass(float* M);
   virtual void ComputeElementForces(const float* U, float* F);
   virtual void ComputeElementPressure(const float* U, float* Pa){;} // Not used in this element
   virtual void ComputeModifiedElementForces(const float* Pa, const float* Va, float* F){;} // Not used in this element
   virtual void GetDhDx(float* Dh);
   virtual void GetStress(float* S);
   virtual void GetStrain(float* E);
   virtual void GetVMStress(float* SVM);
   virtual void GetVMStrain(float* EVM);
   virtual void ComputeStress(float* U);
   virtual float GetStrainEnergy(void);
   virtual void ComputeModifiedStress(float* Pa, float* Va){;} // Not used in this element
   virtual void SetMaterialParams(std::vector<float> params);
   virtual void SetGeometry(tledMesh* Mesh);

private:
   void FillSv();
   void ComputeDhDx();
   
   std::vector<int> EInd;	// Element node indices
   float DhDx[4][3];	// Shape function global derivatives
   float u[4][3];		// Element nodal displacements
   float SPK[3][3];	// 2nd Piola-Kirchoff stress
   float Sv[6];		// Vector version of SPK
   float X[3][3];		// Deformation gradient
   float Fe[3][4];		// Element nodal force contributions
   float x[4][3];    // Element nodal coords
};

#endif // tledElementT4_H
