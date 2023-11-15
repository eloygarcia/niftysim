// =========================================================================
// File:       testMovingRigidContactSurface.cpp
// Purpose:    Moving rigid contact surface unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledUnitTest.h"
#include "tledModel.h"
#include "tledMovingRigidContactSurfaceCPU.h"
#include "tledUnstructuredContactManager.h"
#include "tledSolverCPU.h"

#include <cstdlib>
#include <string>

static void _TestTranslation(const std::string &xmlPath, const int surfaceIndex) {
  using namespace tledVectorArithmetic;

  const int historySize = 20;
  const int totalNumSteps = 1000;

  tledModel model(const_cast<char*>(xmlPath.c_str()));
  tledSolverCPU solver;

  solver.Init(&model);

  {
    tledUnstructuredContactManager surfMan(model, solver, false);
    tledMovingRigidContactSurfaceT3CPU &r_surf = surfMan.GetRigidSurface<tledMovingRigidContactSurfaceT3CPU>(surfaceIndex);
    float trans[3];
  
    std::fill(trans, trans + 3, 0.f);
    r_surf.SetRotations(trans, trans);

    for (float *p_t = trans; p_t < trans + 3; p_t++) *p_t = 10*(drand48() - drand48());
    r_surf.SetTranslation(trans);

    r_surf.SetHistoryLength(historySize);
    r_surf.SetTotalNumberOfSteps(totalNumSteps);

    for (int s = 0; s < totalNumSteps; s++) {
      tledUnitTestAssert(r_surf.GetCurrentStep() == s);

      if (s%5 == 0) {
	for (int t = 0; t < 100; t++) {
	  const int nInd = std::rand()%r_surf.GetNumberOfNodes();
	  const float *x0 = &r_surf.GetAllNodeCoordinates0()[3*nInd];
	  const float *xO = &r_surf.GetAllOldNodeCoordinates()[3*nInd];
	  const float *xN = &r_surf.GetAllNodeCoordinates()[3*nInd];

	  float dx;

	  /* Check node positions */
	  for (int i = 0; i < 3; i++) {
	    if (s < historySize) {
	      tledUnitTestAssert(x0[i] == xO[i]);	    
	      dx = s*trans[i]/totalNumSteps;
	    } else {
	      dx = historySize*trans[i]/totalNumSteps;
	    }

	    tledUnitTestAssert(std::fabs(dx - (xN[i] - xO[i])) <= 1e-4f*std::fabs(dx));	    
	    tledUnitTestAssert(std::fabs(x0[i] + s*trans[i]/totalNumSteps - xN[i]) <= 1e-4f*std::fabs(s*trans[i]/totalNumSteps));	    	  
	  }
	}

	/* Check invariance of normals */
	for (int t = 0; t < 10; t++) {
	  const int fInd = std::rand()%r_surf.GetNumberOfFacets();

	  float diff[3];

	  tledUnitTestAssert(Norm(Sub(diff, r_surf.GetFacetNormal(fInd), r_surf.GetOldFacetNormal(fInd))) < 1e-4f);
	}


	/* Check projection operators */
	for (int t = 0; t < 10; t++) {
	  const int fInd = std::rand()%r_surf.GetNumberOfFacets();
	  const int *facet = r_surf.GetFacet(fInd).NodeIndices;
	  const float deltaN = 0.1;

	  float x[3] = {0.f, 0.f, 0.f}, xi[3] = {-1e10f, -1e10f, -1e10f}, refXi[3] = {0.3333333f, 0.3333333f, 0.f}, diff[3];
	
	  for (int v = 0; v < 3; v++) Add(x, x, r_surf.GetNodeCoordinates(facet[v]));
	  ScalarDiv(x, 3.f);

	  r_surf.ProjectOntoFacet(xi, x, fInd);
	  tledUnitTestAssert(Norm(Sub(diff, xi, refXi)) < 1e-4f);

	  refXi[2] = deltaN;
	  Add(x, x, ScalarMul(diff, r_surf.GetFacetNormal(fInd), deltaN));
	  r_surf.ProjectOntoFacet(xi, x, fInd);
	  tledUnitTestAssert(Norm(Sub(diff, xi, refXi)) < 1e-4f);
	}
      }

      r_surf.Update();
    }
  }
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestTranslation(tledUnitTest::GetResourcePath("moving_wall.xml"), 0);

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
