// =========================================================================
// File:       testSolver.cpp
// Purpose:    tledSolver unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledUnitTest.h"
#include "tledSolverCPU.h"
#include "tledContactManager.h"
#include "tledMesh.h"
#include "tledVectorArithmetic.h"
#include "tledMSHMeshLoader.h"
#ifdef _GPU_
#include "tledCUDAUnitTest.h"
#include "tledSolverGPU.h"
#include "tledSolverGPU_ROM.h"
#endif

#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include <iterator>

using namespace std;

static void _TestAccumulation(const vector<int> inds[], const vector<float> mags[], const float *solverVals, const int numNodes) {
  vector<float> refVals(3*numNodes, 0.0f);

  for (int c = 0; c < 3; c++) for (int n = 0; n < (int)inds[c].size(); n++) {
      refVals[3*inds[c][n]+c] += mags[c][n];
    }

  for (int i = 0; i < 3*numNodes; i++) {
    tledUnitTestAssert(fabs(solverVals[i] - refVals[i]) < 1e-5f);
  }
}

static void _GenRandomIndsAndMags(vector<int> *pv_inds, vector<float> *pv_mags, const int numNodes) {
  for (int c = 0; c < 3; c++) {
    int numC = rand()%(2*numNodes);

    pv_inds[c].resize(numC);
    for (vector<int>::iterator i_i = pv_inds[c].begin(); i_i < pv_inds[c].end(); i_i++) *i_i = rand()%numNodes;
    pv_mags[c].resize(pv_inds[c].size());
    for (vector<float>::iterator i_m = pv_mags[c].begin(); i_m < pv_mags[c].end(); i_m++) *i_m = 10*(float)(drand48() - drand48());
  }
}

static void _TestForces(tledSolver &r_solver) {
  const int numNodes = r_solver.GetMesh()->GetNumNodes();

  vector<int> inds[3], emptyI;
  vector<float> mags[3], solverR(3*numNodes, numeric_limits<float>::quiet_NaN()), emptyM;  

  _GenRandomIndsAndMags(inds, mags, numNodes);
  r_solver.SetDisps(&emptyI, &emptyM, &emptyI, &emptyM, &emptyI, &emptyM);  
  r_solver.SetExtForces(&inds[0], &mags[0], &inds[1], &mags[1], &inds[2], &mags[2]);  
  r_solver.PerformStep();
  r_solver.GetAllExtForces(&solverR.front());
  _TestAccumulation(inds, mags, &solverR.front(), numNodes);
}

static void _TestDisps(tledSolver &r_solver) {
  const int numNodes = r_solver.GetMesh()->GetNumNodes();

  vector<int> inds[3], emptyI;
  vector<float> mags[3], emptyM;

  _GenRandomIndsAndMags(inds, mags, numNodes);
  r_solver.SetExtForces(&emptyI, &emptyM, &emptyI, &emptyM, &emptyI, &emptyM);  
  r_solver.SetDisps(&inds[0], &mags[0], &inds[1], &mags[1], &inds[2], &mags[2]);  
  r_solver.PerformStep();  
  r_solver.PrepareOutput();
  _TestAccumulation(inds, mags, r_solver.GetAllDisps(), numNodes);
}

template <class TSolver>
static void _TestBC(const string &xmlPath) {
  tledModel model(xmlPath.c_str());
  tledContactManager cMgr(&model);
  TSolver solver;
  vector<int> empty;

  solver.Init(&model);
  solver.SetContactManager(&cMgr);
  solver.SetFixed(&empty, &empty, &empty);
  _TestForces(solver);
  _TestDisps(solver);
}

static float _ComputeElementMass(const float x[], const std::vector<int> &vtcs, const float rho) {
  using namespace tledVectorArithmetic;

  float m;

  if (vtcs.size() == 4) {
    float J[3][3], v;

    for (int e = 0; e < 3; e++) Sub(J[e], x + 3*vtcs[e+1], x + 3*vtcs[0]);
    MatDet33(J, &v);
    m = rho*v/6;
  } else {
    float v = 0.f, e[3];

    /* assumes cube/rectangular elements */
    v = Norm(Sub(e, x + 3*vtcs[1], x + 3*vtcs[0]));
    v *= Norm(Sub(e, x + 3*vtcs[3], x + 3*vtcs[0]));
    v *= Norm(Sub(e, x + 3*vtcs[4], x + 3*vtcs[0]));
    m = v*rho;
  }

  return m;
}

template <class TSolver>
static void _TestInhomogeneousMassRectangularStructuredBars(const std::string &meshPath, const char meshType[]) {
  tledUnitTestXMLWriter xmlWriter;
  float refMass = 0.f;

  {
    static const float sectorDensities[3] = {0.5f, 1.75f, 3.25f};
    static const float splitRatios[2] = {0.5f, 0.5f};

    tledMesh mesh = tledUnitTest::LoadMSHMesh(meshPath, meshType);
    std::vector<int> elSets[3];
    float bounds[3][2];
    bool wasSuccess = true;
    
    for (int c = 0; c < 3; c++) bounds[c][1] = -(bounds[c][0] = std::numeric_limits<float>::max());
    for (float const *pc_x = mesh.GetAllNodeCds(); pc_x < mesh.GetAllNodeCds() + 3*mesh.GetNumNodes();) {
      for (int c = 0; c < 3; c++, pc_x++) {
	bounds[c][0] = std::min(*pc_x, bounds[c][0]);
	bounds[c][1] = std::max(*pc_x, bounds[c][1]);
      }
    }
  
    for (int e = 0; e < mesh.GetNumEls(); e++) {
      std::vector<int> vtcs = mesh.GetElNodeInds(e);
      std::vector<int>::const_iterator ic_n;

      assert(vtcs.size() == 4 || vtcs.size() == 8);      
      for (ic_n = vtcs.begin(); ic_n < vtcs.end(); ic_n++) {
	if (mesh.GetAllNodeCds()[*ic_n*3] > (splitRatios[0] + 1e-3f)*bounds[0][1] + (1 - splitRatios[0])*bounds[0][0]) break;
      }

      if (ic_n < vtcs.end()) {
	for (ic_n = vtcs.begin(); ic_n < vtcs.end(); ic_n++) {
	  if (mesh.GetAllNodeCds()[*ic_n*3+1] > (1e-3f + splitRatios[1])*bounds[1][1] + (1 - splitRatios[1])*bounds[1][0]) break;
	}	

	if (ic_n < vtcs.end()) {
	  elSets[2].push_back(e);
	} else {
	  elSets[1].push_back(e);	
	}
      } else {
	elSets[0].push_back(e);	
      }
    }

    refMass = 0.f;
    for (int es = 0; es < 3; es++) {
      for (std::vector<int>::const_iterator ic_e = elSets[es].begin(); ic_e < elSets[es].end(); ic_e++) {
	refMass += _ComputeElementMass(mesh.GetAllNodeCds(), mesh.GetElNodeInds(*ic_e), sectorDensities[es]);
      }
    }
    assert(int(elSets[0].size() + elSets[1].size() + elSets[2].size()) == mesh.GetNumEls());

    wasSuccess &= xmlWriter.StartXML();
    xmlWriter.GetFileStream() << "\t<MSHMesh Type=\"" << meshType << "\">\n"
			      << "\t\t" << meshPath << std::endl
			      << "\t</MSHMesh>\n";
    for (int es = 0; es < 3; es++) {
      assert(elSets[es].size() > 0);
      xmlWriter.GetFileStream() << "\t<ElementSet Size=\"" << elSets[es].size() << "\">\n"
				<< "\t\t<Material Type=\"NH\">\n"      
				<< "\t\t\t<ElasticParams NumParams=\"2\">100 1000</ElasticParams>\n"
				<< "\t\t\t<Density>" << sectorDensities[es] << "</Density>\n"
				<< "\t\t</Material>\n";
      std::copy(elSets[es].begin(), elSets[es].end(), std::ostream_iterator<int>(xmlWriter.GetFileStream(), "\n"));
      xmlWriter.GetFileStream() << "\t</ElementSet>\n";
    }

    xmlWriter.GetFileStream() << "<SystemParams>\n"
			      << "\t<TimeStep>1e-3f</TimeStep>\n"
			      << "\t<TotalTime>2</TotalTime>\n"
			      << "\t<DampingCoeff>3</DampingCoeff>\n"
			      << "\t<HGKappa>0.075</HGKappa>\n"
			      << "</SystemParams>\n";

    wasSuccess &= !xmlWriter.GetFileStream().fail();
    wasSuccess &= xmlWriter.CloseXML();
    assert(wasSuccess);
  }

  {
    tledModel model(const_cast<char*>(xmlWriter.GetFilePath().c_str()));
    TSolver solver;
    std::vector<float> solverMass(model.GetMesh()->GetNumNodes());
    float testMass = 0.f;

    solver.Init(&model);
    solver.GetMassVector(&solverMass.front());
    for (std::vector<float>::const_iterator ic_m = solverMass.begin(); ic_m < solverMass.end(); ic_m++) testMass += *ic_m;

    tledUnitTestAssert(std::fabs(testMass - refMass) < 1e-2f*refMass);
  } 
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestBC<tledSolverCPU>(tledUnitTest::GetResourcePath("moving_wall.xml"));
  _TestInhomogeneousMassRectangularStructuredBars<tledSolverCPU>(tledUnitTest::GetMeshPath("box8.msh"), "H8");
  _TestInhomogeneousMassRectangularStructuredBars<tledSolverCPU>(tledUnitTest::GetMeshPath("box.msh"), "T4");

#ifdef _GPU_
  tledCUDAUnitTest::InitCUDATests();
  _TestBC<tledSolverGPU>(tledUnitTest::GetResourcePath("moving_wall.xml"));
  _TestInhomogeneousMassRectangularStructuredBars<tledSolverGPU>(tledUnitTest::GetMeshPath("box8.msh"), "H8");
  _TestInhomogeneousMassRectangularStructuredBars<tledSolverGPU>(tledUnitTest::GetMeshPath("box.msh"), "T4");
  tledCUDAUnitTest::FinishCUDATests();
#endif

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
