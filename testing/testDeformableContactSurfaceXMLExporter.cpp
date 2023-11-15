// =========================================================================
// File:       testDeformableContactSurfaceXMLExporter.cpp
// Purpose:    Unit test for VTK mesh loader
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledUnitTest.h"
#include "tledDeformableContactSurfaceXMLExporter.h"
#include "tledModel.h"
#include "xmlParser.h"
#include "tledSolverCPU.h"
#include "tledUnstructuredContactManager.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"

#ifdef GPU_GP_CONTACT
#include "tledCUDAUnitTest.h"
#include "tledSolverGPU.h"
#include "tledDeformableContactSurfaceGPU.h"
#endif

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>

template <class TSurface>
void _TestXMLNodes(TSurface &r_surface) {
  tledDeformableContactSurfaceXMLExporter<TSurface> exporter;
  XMLNode root, testChild;

  exporter.SetInput(r_surface);
  exporter.Export();

  root = exporter.GetRootNode();

  tledUnitTestAssert(root.nChildNode("Nodes") == 1);
  tledUnitTestAssert(root.nChildNode("Elements") == 1);
  tledUnitTestAssert(root.nChildNode("Edges") == 1);
  tledUnitTestAssert(root.nChildNode("FacetEdgeIndices") == 1);
  tledUnitTestAssert(root.nChildNode("Volume2SurfaceNodeMap") == 1);
  tledUnitTestAssert(root.nChildNode("Surface2VolumeNodeMap") == 1);
  tledUnitTestAssert(root.nChildNode("NumberOfNodeFacetIndices") == 1);
  tledUnitTestAssert(root.nChildNode("NodeFacetIndices") == 1);

  testChild = root.getChildNode("Nodes");
  
  {
    std::vector<float> xmlN = GetXMLTextAsVector<float>(testChild);
    float const *pc_n = r_surface.GetAllNodeCoordinates();
    std::ostringstream oss;
    
    tledUnitTestAssert((int)xmlN.size() == 3*r_surface.GetNumberOfNodes());

    oss << r_surface.GetNumberOfNodes();
    tledUnitTestAssert(testChild.getAttribute("NumNodes") == oss.str());

    oss.str("");
    oss << 3;
    tledUnitTestAssert(testChild.getAttribute("DOF") == oss.str());

    for (std::vector<float>::const_iterator ic_n = xmlN.begin(); ic_n < xmlN.end(); ic_n++, pc_n++) {
      tledUnitTestAssert(std::fabs(*ic_n - *pc_n) <= std::fabs(*pc_n)*1e-3);
    }
  }

  testChild = root.getChildNode("Elements");
  {
    std::vector<int> xmlF = GetXMLTextAsVector<int>(testChild);
    typename std::vector<typename TSurface::Facet>::const_iterator ic_f = r_surface.GetAllFacets().begin();
    std::ostringstream oss;
    
    tledUnitTestAssert(testChild.getAttribute("Type") == (TSurface::Facet::NumberOfVertices == 3? std::string("T3") : std::string("Q4")));

    oss << r_surface.GetNumberOfFacets();
    tledUnitTestAssert(oss.str() == testChild.getAttribute("NumEls"));

    tledUnitTestAssert((int)xmlF.size() == TSurface::Facet::NumberOfVertices*r_surface.GetNumberOfFacets());
    for (std::vector<int>::const_iterator ic_ni = xmlF.begin(); ic_ni < xmlF.end(); ic_ni += TSurface::Facet::NumberOfVertices, ic_f++) {
      tledUnitTestAssert(std::equal(ic_ni, ic_ni + TSurface::Facet::NumberOfVertices, ic_f->NodeIndices));
    }
  }

  testChild = root.getChildNode("Edges");
  {
    std::vector<int> xmlE = GetXMLTextAsVector<int>(testChild);
    std::vector<std::pair<int, int> >::const_iterator ic_e = r_surface.GetAllEdges().begin();
    
    tledUnitTestAssert((int)xmlE.size() == 2*r_surface.GetNumberOfEdges());
    for (std::vector<int>::const_iterator ic_te = xmlE.begin(); ic_te < xmlE.end(); ic_te += 2, ic_e++) {
      tledUnitTestAssert(*ic_te == ic_e->first && *(ic_te + 1) == ic_e->second);
    }
  }  

  testChild = root.getChildNode("NumberOfNodeFacetIndices");
  {
    std::vector<int> xmlB = GetXMLTextAsVector<int>(testChild);
    
    tledUnitTestAssert((int)xmlB.size() == r_surface.GetNumberOfNodes());
    for (int n = 0; n < r_surface.GetNumberOfNodes(); n++) {
      tledUnitTestAssert(xmlB[n] == r_surface.GetNumberOfNodeFacets(n));
    }
  }

  testChild = root.getChildNode("NodeFacetIndices");
  {
    std::vector<int> xmlNF = GetXMLTextAsVector<int>(testChild);
    std::vector<int>::const_iterator ic_nf = xmlNF.begin();
    int tn = 0;

    for (int n = 0; n < r_surface.GetNumberOfNodes(); n++) tn += r_surface.GetNumberOfNodeFacets(n);
    tledUnitTestAssert((int)xmlNF.size() == tn);

    for (int n = 0; n < r_surface.GetNumberOfNodes(); n++) {      
      tledUnitTestAssert(std::equal(ic_nf, ic_nf + r_surface.GetNumberOfNodeFacets(n), r_surface.GetNodeFacetIndices(n)));
      ic_nf += r_surface.GetNumberOfNodeFacets(n);
    }
  }

  root.deleteNodeContent();
}

template <class TSolver, class TSurfaceAPI, class TSurface, class TMembraneSurface, const bool t_useGPU>
static void _TestXMLNodes(const std::string &inSimFile) {
  tledModel model(inSimFile.c_str());
  TSolver solver;

  solver.Init(&model);

  {
    tledUnstructuredContactManager surfaces(model, solver, t_useGPU);
    TSurfaceAPI &r_surface = surfaces.GetDeformableSurface<TSurfaceAPI>();

    assert(&r_surface != NULL);
    if (r_surface.GetNumberOfFacetVertices() == 3 && (model.GetNumberOfShellElementSets() == 0 || model.DoShellUseMeshSurface())) {
      _TestXMLNodes(dynamic_cast<TSurface&>(r_surface));
    } else if (r_surface.GetNumberOfFacetVertices() == 3 && model.GetNumberOfShellElementSets() > 0 && !model.DoShellUseMeshSurface()) {
#ifdef GPU_GP_CONTACT
      if (t_useGPU) {
	tledFatalNotYetImplementedError;
      }
#endif
      _TestXMLNodes(dynamic_cast<TMembraneSurface&>(r_surface));
    } else {
      tledFatalError("Incompatible contact surface.");
    }
  }
}

static void _TestXMLNodesCPU(const std::string &inSimFile) {
  _TestXMLNodes<tledSolverCPU, tledDeformableContactSurfaceCPU, tledDeformableContactSurfaceT3CPU, tledDeformableMembraneContactSurfaceT3CPU, false>(inSimFile);
}

#ifdef GPU_GP_CONTACT
static void _TestXMLNodesGPU(const std::string &inSimFile) {
  _TestXMLNodes<tledSolverGPU, tledDeformableContactSurfaceGPU, tledDeformableContactSurfaceT3GPU, tledDeformableContactSurfaceT3GPU, true>(inSimFile);
}
#endif

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestXMLNodesCPU(tledUnitTest::GetMeshPath("ball_ball.xml"));
  _TestXMLNodesCPU(tledUnitTest::GetMeshPath("two_cubes_sheet.xml"));

#ifdef GPU_GP_CONTACT
  tledCUDAUnitTest::InitCUDATests();

  _TestXMLNodesGPU(tledUnitTest::GetMeshPath("ball_ball.xml"));

  tledCUDAUnitTest::FinishCUDATests();
#endif

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
