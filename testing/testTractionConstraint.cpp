// =========================================================================
// File:       testTractionConstraint.cpp
// Purpose:    TractionConstraint unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledUnitTest.h"
#include "tledSimulator.h" /* Do NOT include tledTractionConstraint, this is a test in and of itself! */
#include "tledMeshSurface.h"

#include <cmath>
#include <limits>

template <const int t_numFacetVtcs>
static void _ExtractCubeSideFaces(std::vector<int> &r_faces, float &r_area, const tledMesh &mesh, const float cVal, const int axis) {
  const float *x = mesh.GetAllNodeCds();
  const tledMeshSurface<t_numFacetVtcs> surface(mesh);

  r_faces.clear();
  r_area = 0.f;
  for (int f = 0; f < surface.GetNumberOfFacets(); f++) {
    int v;

    for (v = 0; v < t_numFacetVtcs && std::fabs(x[3*surface.GetFacet(f).NodeIndices[v]+axis] - cVal) <= 1e-2f*std::fabs(cVal); v++);
    if (v == t_numFacetVtcs) {
      r_faces.insert(r_faces.end(), surface.GetFacet(f).NodeIndices, surface.GetFacet(f).NodeIndices + t_numFacetVtcs);
      r_area += surface.ComputeFacetArea(f);
    }
  }
}

template <const int t_numFacetVtcs>
static void _TestUnitCubeNoModel() {
  const std::string meshtype = t_numFacetVtcs == 3? "T4" : "H8";
  const std::string meshPath = tledUnitTest::GetMeshPath(t_numFacetVtcs == 3? "box.msh" : "box8.msh");
  const tledMesh cube = tledUnitTest::LoadMSHMesh(meshPath, meshtype.c_str());
  const float *x = cube.GetAllNodeCds();
  
  tledModel *p_model = NULL;
  tledSimulator *p_sim = NULL;

  {
    tledUnitTestXMLWriter xmlWriter;

    xmlWriter.StartXML();
    xmlWriter.GetFileStream() << "\t<MSHMesh Type=\"" << meshtype << "\">\n"
			      << "\t\t" << meshPath << endl
			      << "\t</MSHMesh>\n";
    xmlWriter.GetFileStream() << "\t<ElementSet Size=\"all\">\n"
			      << "\t\t<Material Type=\"NH\">\n"      
			      << "\t\t\t<ElasticParams NumParams=\"2\">100 1000</ElasticParams>\n"
			      << "\t\t\t<Density>10</Density>\n"
			      << "\t\t</Material>\n"
			      << "\t</ElementSet>\n";
    xmlWriter.GetFileStream() << "<SystemParams>\n"
			      << "\t<TimeStep>1e-3f</TimeStep>\n"
			      << "\t<TotalTime>0.01</TotalTime>\n"
			      << "\t<DampingCoeff>1</DampingCoeff>\n"
			      << "\t<HGKappa>0.075</HGKappa>\n"
			      << "</SystemParams>\n";
    if (xmlWriter.GetFileStream().fail() || !xmlWriter.CloseXML()) {
      tledFatalError("Error writing unit test XML");
    }    

    p_model = new tledModel(xmlWriter.GetFilePath().c_str());
    p_sim = new tledSimulator(p_model);
    p_sim->Simulate();
    p_sim->GetSolver()->PrepareOutput();
  }  

  for (int a = 0; a < 3; a++) for (int s = 0; s < 2; s++) {
      const float constTrac[] = {0.f, 0.f, 1.23f};

      float minAVal, maxAVal, area;
      std::vector<float> tractions;
      std::vector<int> faces, nodes;
      
      maxAVal = -(minAVal = std::numeric_limits<float>::max());
      for (float const *pc_x = x + a; pc_x < x + 3*cube.GetNumNodes(); pc_x += 3) {
	maxAVal = std::max(maxAVal, *pc_x);
	minAVal = std::min(minAVal, *pc_x);
      }
      
      if (s == 0) _ExtractCubeSideFaces<t_numFacetVtcs>(faces, area, cube, minAVal, a);
      else _ExtractCubeSideFaces<t_numFacetVtcs>(faces, area, cube, maxAVal, a);
      
      tractions.reserve(int(faces.size()/t_numFacetVtcs)*3);
      for (int f = 0; f < int(faces.size()/t_numFacetVtcs); f++) tractions.insert(tractions.end(), constTrac, constTrac + 3);
      
      {
	tledTractionConstraint constraint(*p_sim->GetSolver(), t_numFacetVtcs == 4? 0 : 1, faces, tractions, RAMP);

	nodes = tledHelper::MakeSortedUnique(faces);
	tledUnitTestAssert(std::equal(nodes.begin(), nodes.end(), constraint.GetForceInd(a)->begin()));

	std::fill(p_sim->GetSolver()->GetAllDisps(), p_sim->GetSolver()->GetAllDisps() + p_sim->GetModel()->GetMesh()->GetNumNodes()*3, 0.f);
	for (int t = 0; t < 2; t++) {
	  const float refVal = area*constTrac[2]*float(t + 1)/2;

	  std::vector<float> *p_f;
	  float tSum = 0.f;
	
	  p_f = constraint.GetForceVal(2, t, 0.1, 2*0.1);
	  for (std::vector<float>::const_iterator ic_f = p_f->begin(); ic_f < p_f->end(); ic_f++) tSum += *ic_f;
	  tledUnitTestAssert((refVal + tSum) <= std::max(1e-2f*std::fabs(refVal), 1e-4f));
	}
      }
    }
    
  delete p_sim;
  delete p_model;
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestUnitCubeNoModel<3>();
  _TestUnitCubeNoModel<4>();

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
