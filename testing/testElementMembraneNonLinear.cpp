// =========================================================================
// File:       testElementMembraneNonLinear.cpp
// Purpose:    Membrane element unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledUnitTest.h"
#include "tledCUDAUnitTest.h"
#include "tledShellPatchTest.h"
#include "tledCUDAHelpers.h"

int main(void) {
  using namespace tledShellPatchTest;

  tledUnitTest::InitUnitTest();

  TestForceStretch(tledUnitTest::GetResourcePath("membrane_force_stretch_nonlinear.xml"), false, 1.f, 1.f, 1.f, 1e3f, 25.f);
  TestMass(tledUnitTest::GetResourcePath("membrane_static_nonlinear.xml"), 1.f);
  TestRigid(tledUnitTest::GetResourcePath("membrane_static_nonlinear.xml"), 0.1f);
  TestStretch(tledUnitTest::GetResourcePath("membrane_stretch_nonlinear.xml"), false, 0.01f, 1.f, 0.5f);
  TestStretch(tledUnitTest::GetResourcePath("membrane_stretch_nonlinear_irregular.xml"), false, 0.01f, 1.f, 0.5f);
  TestShear(tledUnitTest::GetResourcePath("membrane_shear_nonlinear.xml"), false, 0.05f);

#ifdef _GPU_
  tledCUDAUnitTest::InitCUDATests();

  TestForceStretch(tledUnitTest::GetResourcePath("membrane_force_stretch_nonlinear.xml"), true, 1.f, 1.f, 1.f, 1e3f, 25);
  TestStretch(tledUnitTest::GetResourcePath("membrane_stretch_nonlinear_irregular.xml"), true, 0.01f, 1.f, 0.5f);
  TestShear(tledUnitTest::GetResourcePath("membrane_shear_nonlinear.xml"), true, 0.05f);

  tledCUDAUnitTest::FinishCUDATests();
#endif

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
}
