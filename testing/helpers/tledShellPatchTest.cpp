// =========================================================================
// File:       tledShellPatchTest.cpp
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

#ifdef _GPU_
#undef _GPU_
#endif

#ifdef GPU_GP_CONTACT
#undef GPU_GP_CONTACT
#endif

#include "tledShellPatchTest.h"
#include "tledVectorArithmetic.h"
#include "tledModel.h"
#include "tledSimulator.h"
#include "tledShellSolverCPU.h"
#include "tledUnitTest.h"

#include <vector>

namespace tledShellPatchTest {
  using namespace tledVectorArithmetic;

  void TestRigid(const std::string &xmlPath, const float t0) {
    tledModel model(const_cast<char*>(xmlPath.c_str()));
    tledSimulator sim(&model, false, false, false);
    tledShellSolver &r_shellSolver = sim.GetSolver()->GetShellSolver();
    std::vector<float> U, F;
  
    sim.Simulate();
    sim.GetSolver()->PrepareOutput();
    U.resize(3*r_shellSolver.GetSurface().GetNumberOfNodes());
    F.resize(3*r_shellSolver.GetSurface().GetNumberOfNodes());

    for (int t = 0; t < 2; t++) {
      switch (t) {
      case 0: {
	const float dx = 0.3f, dy = 0.7f, dz = 0.8f;

	for (std::vector<float>::iterator i_u = U.begin(); i_u < U.end();) {
	  *(i_u++) = dx;
	  *(i_u++) = dy;
	  *(i_u++) = dz;
	}

	break;
      }

      case 1: {
	const float cor[] = {0.5f, 0.5f, 0.f};
	const float rot[] = {0.3f, 0.2f, 0.1f};

	for (int n = 0; n < r_shellSolver.GetSurface().GetNumberOfNodes(); n++) {
	  const float *x0 = r_shellSolver.GetSurface().GetNodeCoordinates(n);

	  float xCnt[3], xNew[3];
	  float nOld, nNew;

	  Sub(xCnt, x0, cor);
	  nOld = Norm(xCnt);

	  xNew[0] = std::cos(rot[0])*xCnt[0] - std::sin(rot[0])*xCnt[1];
	  xNew[1] = std::sin(rot[0])*xCnt[0] + std::cos(rot[0])*xCnt[1];
	  xNew[2] = xCnt[2];
	
	  xCnt[0] = std::cos(rot[1])*xNew[0] - std::sin(rot[1])*xNew[2];
	  xCnt[1] = xNew[1];
	  xCnt[2] = std::sin(rot[1])*xNew[0] + std::cos(rot[1])*xNew[2];

	  xNew[0] = xCnt[0];	
	  xNew[1] = std::cos(rot[0])*xCnt[1] - std::sin(rot[0])*xCnt[2];
	  xNew[2] = std::sin(rot[0])*xCnt[1] + std::cos(rot[0])*xCnt[2];

	  nNew = Norm(xNew);
	  assert(std::fabs(nOld - nNew) < 1e-3f);

	  Sub(&U.front() + 3*n, Add(&U.front() + 3*n, xNew, cor), x0);
	}

	break;
      }
      
      }

      {
	float avgF = 0;
      
	std::fill(F.begin(), F.end(), 0.0f);
	static_cast<tledShellSolverCPU&>(r_shellSolver).ComputeNewForces(&F.front(), &U.front());
	for (std::vector<float>::const_iterator ic_f = F.begin(); ic_f < F.end(); ic_f++) avgF += *ic_f;

	tledUnitTestAssert(std::fabs(avgF) < 1e-3f);
      }
      
      {
	float avgT = 0.f, stdT = 0.f;
	std::vector<float> elementTs(r_shellSolver.GetNumberOfElementSetElements(0));
    
	r_shellSolver.ComputeElementThicknesses(&elementTs.front(), 0, &U.front());
	for (std::vector<float>::const_iterator ic_t = elementTs.begin(); ic_t < elementTs.end(); ic_t++) avgT += *ic_t;
	avgT /= elementTs.size();

	for (std::vector<float>::const_iterator ic_t = elementTs.begin(); ic_t < elementTs.end(); ic_t++) stdT += (*ic_t - avgT)*(*ic_t - avgT);
	stdT = std::sqrt(stdT/elementTs.size());

	tledUnitTestAssert(std::fabs(avgT - t0) < 1e-3f);
	tledUnitTestAssert(stdT < 1e-3f);
      }
    }
  }

  void TestShear(const std::string &path, const bool useGPU, const float deltaY) {
    tledModel model(const_cast<char*>(path.c_str()));
    tledSimulator sim(&model, useGPU, false, false);
    tledShellSolver &r_shellSolver = sim.GetSolver()->GetShellSolver();
    float const *pc_U = NULL;
  
    sim.Simulate();

    sim.GetSolver()->PrepareOutput();
    pc_U = sim.GetSolver()->GetAllDisps();

    {
      const tledSurface &mesh = r_shellSolver.GetSurface();

      for (int nInd = 0; nInd < mesh.GetNumberOfNodes(); nInd++) {
	const float x = mesh.GetNodeCoordinates(nInd)[0];
	const float dy = deltaY*x;

	if (mesh.GetNodeCoordinates(nInd)[1] == 1 || mesh.GetNodeCoordinates(nInd)[1] == 0) {
	  tledUnitTestAssert(std::fabs(*(pc_U + 3*nInd + 1) - dy) < 1e-3f);
	}
      }
    }
  }

  void TestStretch(const std::string &xmlPath, const bool useGPU, const float deltaL, const float t0, const float nu) {
    tledModel model(const_cast<char*>(xmlPath.c_str()));
    tledSimulator sim(&model, useGPU, false, false);
    tledShellSolver &r_shellSolver = sim.GetSolver()->GetShellSolver();    
    float const *pc_U = NULL;
  
    sim.Simulate();

    sim.GetSolver()->PrepareOutput();
    pc_U = sim.GetSolver()->GetAllDisps();

    {
      const tledSurface &mesh = r_shellSolver.GetSurface();

      float deltaYTop = 0, deltaYBottom = 0;
      int numTop = 0, numBot = 0;

      for (int nInd = 0; nInd < mesh.GetNumberOfNodes(); nInd++) {
	if (mesh.GetNodeCoordinates(nInd)[1] == 1) deltaYTop -= *(pc_U + 3*nInd + 1), numTop += 1;
	if (mesh.GetNodeCoordinates(nInd)[1] == 0) deltaYBottom += *(pc_U + 3*nInd + 1), numBot += 1;
      }

      deltaYTop /= numTop;
      deltaYBottom /= numBot;

      assert(deltaYBottom >= 0 && deltaYTop >= 0);
    
      tledUnitTestAssert(std::fabs(deltaYTop - nu*0.5f*deltaL) < 1e-3f);
      tledUnitTestAssert(std::fabs(deltaYBottom - nu*0.5f*deltaL) < 1e-3f);
    }
  
    {
      std::vector<float> elementTs(r_shellSolver.GetNumberOfElementSetElements(0));
    
      r_shellSolver.ComputeElementThicknesses(&elementTs.front(), 0, pc_U);
      //for (std::vector<float>::const_iterator ic_t = elementTs.begin(); ic_t < elementTs.end(); ic_t++) tledUnitTestAssert(std::fabs(*ic_t - t0*(1 - nu*deltaL)) < 1e-3);
    }
  }

  void TestCooksMembrane(const std::string &xmlPath, const bool useGPU) {
    /* Reference value from Liu et al.: A novel alpha finite element method (aFEM) for exact solution to mechanics problems using triangular and tetrahedral elements */
    static const float refCentralDisplacement = 23.9642f;

    tledModel model(const_cast<char*>(xmlPath.c_str()));
    tledSimulator sim(&model, useGPU, false, false);
    tledShellSolver &r_shellSolver = sim.GetSolver()->GetShellSolver();    
    float const *pc_U = NULL;
    int refCornerInd = -1, sndRefCornerInd = -1;
    float minDist = std::numeric_limits<float>::max(), sndMinDist = std::numeric_limits<float>::max();
  
    sim.Simulate();

    sim.GetSolver()->PrepareOutput();
    pc_U = sim.GetSolver()->GetAllDisps();

    {
      const tledSurface &mesh = r_shellSolver.GetSurface();
      const float cp[] = {48.f, 52.f, 0.f};

      for (int nInd = 0; nInd < mesh.GetNumberOfNodes(); nInd++) {
	const float *x = mesh.GetNodeCoordinates(nInd);

	float d[3], dn;

	if ((dn = Norm(Sub(d, cp, x))) < minDist) {
	  refCornerInd = nInd;
	  minDist = dn;
	} else if (dn < sndMinDist) {
	  sndRefCornerInd = nInd;
	  sndMinDist = dn;
	}
      }
    }

    tledUnitTestAssert(std::fabs(refCentralDisplacement - pc_U[3*refCornerInd+1]*sndMinDist/(minDist + sndMinDist) - pc_U[3*sndRefCornerInd+1]*minDist/(minDist + sndMinDist)) < 2e-2f*refCentralDisplacement);
  }    

  void TestMass(const std::string &xmlPath, const float refMass) {
    tledModel model(const_cast<char*>(xmlPath.c_str()));
    tledSimulator sim(&model, false, false, false);
    tledSolverCPU &r_solver = *static_cast<tledSolverCPU*>(sim.GetSolver());    
    std::vector<float> mass;
    float total = 0;

    mass.resize(r_solver.GetMesh()->GetNumNodes());
    r_solver.GetMassVector(&mass.front());
    for (std::vector<float>::const_iterator ic_m = mass.begin(); ic_m < mass.end(); ic_m++) total += *ic_m;

    tledUnitTestAssert(std::fabs(total - refMass) < 1e-4f*refMass);
  }

  void TestCantileverBeam(const std::string &xmlPath, const bool useGPU, const float L, const float D, const float E, const float nu, const float F) {
    tledModel model(const_cast<char*>(xmlPath.c_str()));
    tledSimulator sim(&model, useGPU, false, false);
    tledShellSolver &r_shellSolver = sim.GetSolver()->GetShellSolver();    
    float const *pc_U = NULL;
  
    sim.Simulate();

    sim.GetSolver()->PrepareOutput();
    pc_U = sim.GetSolver()->GetAllDisps();

    {
      const tledSurface &mesh = r_shellSolver.GetSurface();

      for (int nInd = 0; nInd < mesh.GetNumberOfNodes(); nInd++) if (std::fabs(L - mesh.GetNodeCoordinates(nInd)[0]) < 1e-2f*L) {
	const float uy = -2*F*(1.0f/4.0f*(4 + 5*nu)*D*D*L + 2*L*L*L)/(E*D*D*D);

	tledUnitTestAssert(std::fabs(*(pc_U + 3*nInd + 1) - uy) < 1e-2f*std::fabs(uy));
      }
    }
  }

  void TestForceStretch(const std::string &xmlPath, const bool useGPU, const float L, const float W, const float T, const float E, const float F) {
    const float dL = F*L/(W*E*T);

    tledModel model(const_cast<char*>(xmlPath.c_str()));
    tledSimulator sim(&model, useGPU, false, false);
    tledShellSolver &r_shellSolver = sim.GetSolver()->GetShellSolver();    
    float const *pc_U = NULL;
  
    sim.Simulate();

    sim.GetSolver()->PrepareOutput();
    pc_U = sim.GetSolver()->GetAllDisps();

    {
      const tledSurface &mesh = r_shellSolver.GetSurface();

      for (int nInd = 0; nInd < mesh.GetNumberOfNodes(); nInd++) if (std::fabs(L - mesh.GetNodeCoordinates(nInd)[0]) < 1e-2f*L) {
	tledUnitTestAssert(std::fabs(*(pc_U + 3*nInd) - dL) < 1e-3f*std::fabs(L));
      }
    }
  }
}
