// =========================================================================
// File:       tledElementT4ANP.h
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


#ifndef tledElementT4ANP_H
#define tledElementT4ANP_H

#include "tledElement.h"
#include "tledMesh.h"

class tledElementT4ANP : public tledElement
{
public:
   tledElementT4ANP() {;}
   tledElementT4ANP(tledMesh* Mesh, int ElNum, const char* MatType, float* MatParams, float Density);
   virtual ~tledElementT4ANP();
   
   virtual void ComputeElementMass(float* M);
   virtual void ComputeElementForces(const float* U, float* F){;} // Not used in this element
   virtual void ComputeElementPressure(const float* U, float* Pa);
   virtual void ComputeModifiedElementForces(const float* Pa, const float* Va, float* F);
   virtual void GetDhDx(float* Dh);
   virtual void GetStress(float* S);
   virtual void GetStrain(float* E);
   virtual void GetVMStress(float* SVM);
   virtual void GetVMStrain(float* EVM);
   virtual void ComputeStress(float* U){;} // Not used in this element
   virtual float GetStrainEnergy(void);
   virtual void ComputeModifiedStress(float* Pa, float* Va);
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
   float J;			// Element Jacobian
   float x[4][3];    // Element nodal coords
};

#endif // tledElementT4ANP_H
