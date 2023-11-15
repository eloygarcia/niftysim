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

  TestRigid(tledUnitTest::GetResourcePath("membrane_static_linear.xml"), 0.05f);
  TestStretch(tledUnitTest::GetResourcePath("membrane_stretch_linear.xml"), false, 0.1f, 0.1f, 0.4f);

#ifdef _GPU_
  tledCUDAUnitTest::InitCUDATests();

  TestStretch(tledUnitTest::GetResourcePath("membrane_stretch_linear.xml"), true, 0.1f, 0.1f, 0.4f);

  tledCUDAUnitTest::FinishCUDATests();
#endif

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
}
