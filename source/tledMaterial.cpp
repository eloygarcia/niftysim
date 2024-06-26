// =========================================================================
// File:       tledMaterial.cpp
// Purpose:    Material base class
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


#include "tledMaterial.h"

tledMaterial::tledMaterial()
{
}

tledMaterial::tledMaterial(const float rho) {
  this->SetDensity(rho);
}
