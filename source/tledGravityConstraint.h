// =========================================================================
// File:       tledGravityConstraint.h
// Purpose:    Gravity constraint class. Used for applying gravitational
//             loads.
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    March 2011
// 
// Copyright (c) 2011, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledGravityConstraint_H
#define tledGravityConstraint_H

#include "tledConstraint.h"

/**
 * \brief Constant body-force constraint
 * \ingroup constraints
 */
class tledGravityConstraint : public tledConstraint
{
public:
   tledGravityConstraint();
   tledGravityConstraint(std::vector<int> ind, std::vector<float> mass, float mag, std::vector<float> dir, enum loadShape ls);
   virtual ~tledGravityConstraint();
   
   // Get params
   virtual std::vector<int>& GetInd() {return Ind;}
   virtual const std::vector<int>& GetInd() const {return Ind;}
   float GetAccelMag() {return AccelMag;}
   std::vector<float> GetAccelDir() {return AccelDir;}
   
   // Set params
   void SetInd(std::vector<int> ind) {Ind = ind;}
   void SetAccelMag(float mag) {AccelMag = mag;}
   void SetAccelDir(std::vector<float> dir);
   
   // Methods required by tledConstraintManager
   virtual std::vector<int>* GetForceInd(int dof);
   virtual std::vector<float>* GetForceVal(int dof, int step, double dt, double T);
   
private:
   std::vector<int> Ind;  // Node indices
   std::vector<float> Mass; // Masses of the nodes listed in Ind
   float AccelMag; // Acceleration magnitude
   std::vector<float> AccelDir; // Acceleration direction
   std::vector<float> R;  // Forces for current step
};


#endif // tledGravityConstraint_H
