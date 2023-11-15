// =========================================================================
// File:       testDeformableContactSurfaceCPU.cpp
// Purpose:    tledDeformableContactSurfaceCPU unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    August 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledUnitTest.h"
#include "tledDeformableContactSurfaceCPU.h"

static void _TestNormalsOnSphereCPU() {
  using namespace tledVectorArithmetic;

  const float maxDev = 5.f*tledPi/180.f;

  tledDeformableContactSurfaceT3CPU surface(tledUnitTest::LoadMSHMesh(tledUnitTest::GetMeshPath("sphere.msh"), "T4"));
  float err = 0;
  
  surface.Init();
  for (int nInd = 0; nInd < surface.GetNumberOfNodes(); nInd++) {
    const float *n = surface.GetNodeNormalCached(nInd);
    const float *x = surface.GetNodeCoordinates(nInd);

    float ref[3];

    assert(std::fabs(Norm(x) - 1) <= 1e-2f);
    tledUnitTestAssert(std::fabs(Norm(n) - 1) <= 1e-4f);

    ScalarDiv(ref, x, Norm(x));
    tledUnitTestAssert(std::fabs(ComputeAngleNormalised(n, ref)) < maxDev);
    err += std::fabs(ComputeAngleNormalised(n, ref));
  }

  err /= surface.GetNumberOfNodes();
  std::cout << "Avg. normal error: " << err << std::endl;
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestNormalsOnSphereCPU();

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
