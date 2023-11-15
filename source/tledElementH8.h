// =========================================================================
// File:       tledElementH8.h
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


#ifndef tledElementH8_H
#define tledElementH8_H

#include "tledElement.h"
#include "tledMesh.h"

class tledElementH8 : public tledElement
{
public:
   tledElementH8() {;}
   tledElementH8(tledMesh* Mesh, int ElNum, const char* MatType, float* MatParams, float Density, float kappa);
   virtual ~tledElementH8();
   
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
   void ComputeDhDx_BmatMeth();
   void ComputeVol();
   void ComputeBmat(float B[8][3]);
   void ComputeCIJK(float C[8][8][8]);
   void ComputeHGParams();
   void ComputeHGForces();
   
   int ENum;
   std::vector<int> EInd;
   float detJ;		// Jacobian determinant
   float DhDx[8][3];	// Shape function global derivatives
   float u[8][3];		// Element nodal displacements
   float X[3][3];		// Deformation gradient
   float SPK[3][3];	// 2nd Piola-Kirchoff stress
   float Sv[6];		// Vector version of SPK
   float Fe[3][8];		// Element nodal force contributions
   float x[8][3];    // Element nodal coords
   float HGKappa;    // Hourglass control parameter
   
   float HG[8][8];
   float FHG[8][3];
};

#endif // tledElementH8_H
