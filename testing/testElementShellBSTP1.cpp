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
#include "tledShellSolver.h"
#include "tledModel.h"
#include "tledSimulator.h"
#include "tledElementShellBSTP1.h"
#include "tledCUDAHelpers.h"

int main(void) {
  using namespace tledShellPatchTest;

  tledUnitTest::InitUnitTest();

  TestRigid(tledUnitTest::GetResourcePath("membrane_static_shell.xml"), 0.1f);

  TestShear(tledUnitTest::GetResourcePath("membrane_shear_shell.xml"), false, 0.05f);
  TestForceStretch(tledUnitTest::GetResourcePath("membrane_force_stretch_shell.xml"), false, 1.f, 1.f, 1.f, 1e3f, 25);
  TestCooksMembrane(tledUnitTest::GetResourcePath("membrane_cooks_shell_hires.xml"), false);
  TestMass(tledUnitTest::GetResourcePath("membrane_stretch_shell.xml"), 0.001f);
  TestMass(tledUnitTest::GetResourcePath("membrane_stretch_shell_irregular.xml"), 0.001f);
  TestStretch(tledUnitTest::GetResourcePath("membrane_stretch_shell.xml"), false, 0.05f, 0.1f, 0.3f);
  TestStretch(tledUnitTest::GetResourcePath("membrane_stretch_shell_irregular.xml"), false, 0.05f, 0.1f, 0.4f);

#ifdef _GPU_
  tledCUDAUnitTest::InitCUDATests();

  TestStretch(tledUnitTest::GetResourcePath("membrane_stretch_shell.xml"), true, 0.05f, 0.1f, 0.3f);
  TestForceStretch(tledUnitTest::GetResourcePath("membrane_force_stretch_shell.xml"), true, 1.f, 1.f, 1.f, 1e3f, 25);
  TestShear(tledUnitTest::GetResourcePath("membrane_shear_shell.xml"), true, 0.05f);

  tledCUDAUnitTest::FinishCUDATests();
#endif

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
}
