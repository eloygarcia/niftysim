// =========================================================================
// File:       tledShellPatchTest.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    August 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellPatchTest_H
#define tledShellPatchTest_H

#include <string>

namespace tledShellPatchTest {
  void TestRigid(const std::string &xmlPath, const float t0);
  
  /** Stretching of unit square */
  void TestStretch(const std::string &xmlPath, const bool useGPU, const float deltaL, const float t0, const float nu);

  /** Shearing of unit square */
  void TestShear(const std::string &path, const bool useGPU, const float deltaY);

  void TestForceStretch(const std::string &xmlPath, const bool useGPU, const float L, const float W, const float T, const float E, const float F);

  void TestCooksMembrane(const std::string &xmlPath, const bool useGPU);

  /** Computes the sum of nodal masses */
  void TestMass(const std::string &xmlPath, const float refMass);

  /** Assumes a beam with dimensions [0, L] x [0, D] (centre line at D/2) */
  void TestCantileverBeam(const std::string &xmlPath, const bool useGPU, const float L, const float W, const float E, const float nu, const float F);
}

#endif
