// =========================================================================
// File:       tledElement.h
// Purpose:    Element base class
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


#ifndef tledElement_H
#define tledElement_H

#include "tledMesh.h"
#include "tledMaterial.h"
#include "tledMaterialLE.h"
#include "tledMaterialNH.h"
#include "tledMaterialNHV.h"
#include "tledMaterialTI.h"
#include "tledMaterialTIV.h"
#include "tledMaterialAB.h"
#include "tledMaterialPY.h"

class tledElement
{
public:
   tledElement();
   virtual ~tledElement();
   
   virtual void ComputeElementMass(float* M) = 0;
   virtual void ComputeElementForces(const float* U, float* F) = 0;
   virtual void ComputeElementPressure(const float* U, float* Pa) = 0;
   virtual void ComputeModifiedElementForces(const float* Pa, const float* Va, float* F) = 0;
   virtual void GetDhDx(float* Dh) = 0;	// Get shape function global derivatives
   virtual void GetStress(float* S) = 0;	// Get 2nd Piola-Kirchhoff stress
   virtual void GetStrain(float* E) = 0;	// Get Green-Lagrange strain
   virtual void GetVMStress(float* SVM) = 0;	// Get Von Mises equivalent stress
   virtual void GetVMStrain(float* EVM) = 0;	// Get Von Mises equivalent strain
   virtual void ComputeStress(float* U) = 0;  // Compute the SPK stress for the current configuration
   virtual float GetStrainEnergy(void) = 0;
   virtual void ComputeModifiedStress(float* Pa, float* Va) = 0;
   virtual void SetMaterialParams(std::vector<float> params) = 0;
   virtual void SetGeometry(tledMesh* Mesh) = 0; // For altering the element geometry
   float GetVolume() {return Vol;}
   tledMaterial* GetMaterial() {return Mat;}

protected:
   int NDOF;	// Degrees of freedom per node
   float Vol;	// Element volume
   tledMaterial* Mat;	// Element material
};

#endif // tledElement_H
