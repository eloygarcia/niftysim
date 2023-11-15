// =========================================================================
// File:       testDeformableMembraneContactSurfaceXMLExporter.cpp
// Purpose:    tledDeformableMembraneContactSurfaceXMLExporter unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledUnitTest.h"
#include "tledDeformableMembraneContactSurfaceXMLExporter.h"
#include "tledModel.h"
#include "xmlParser.h"
#include "tledSolverCPU.h"
#include "tledUnstructuredContactManager.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"

/** Only tests membrane specific members */
template <class TSurface>
void _TestXMLNodes(TSurface &r_surface) {
  tledDeformableMembraneContactSurfaceXMLExporter<TSurface> exporter;
  XMLNode root, testChild;
  std::istringstream iss;
  int testVal = -1;

  exporter.SetInput(r_surface);
  exporter.Export();

  root = exporter.GetRootNode();

  tledUnitTestAssert(root.nChildNode("NumberOfMembraneNodes") == 1);
  tledUnitTestAssert(root.nChildNode("NumberOfMembraneFacets") == 1);
  tledUnitTestAssert(root.nChildNode("NumberOfMembraneEdges") == 1);

  tledUnitTestAssert(root.nChildNode("MembraneFacetBaseIndex") == 1);
  tledUnitTestAssert(root.nChildNode("MembraneNodeBaseIndex") == 1);
  tledUnitTestAssert(root.nChildNode("MembraneEdgeBaseIndex") == 1);

  iss.str(root.getChildNode("NumberOfMembraneNodes").getText());
  iss >> testVal;
  tledUnitTestAssert(!iss.fail() && testVal == r_surface.GetNumberOfMembraneNodes());
  iss.clear();

  iss.str(root.getChildNode("NumberOfMembraneFacets").getText());
  iss >> testVal;
  tledUnitTestAssert(!iss.fail() && testVal == r_surface.GetNumberOfMembraneElements());
  iss.clear();

  iss.str(root.getChildNode("NumberOfMembraneEdges").getText());
  iss >> testVal;
  tledUnitTestAssert(!iss.fail() && testVal == r_surface.GetNumberOfMembraneEdges());    
}

static void _TestXMLNodes(const std::string &inSimFile) {
  tledModel model(inSimFile.c_str());
  tledSolverCPU solver;

  solver.Init(&model);

  {
    tledUnstructuredContactManager surfaces(model, solver, false);
    tledDeformableContactSurface &r_surface = surfaces.GetDeformableSurface<tledDeformableContactSurface>();

    assert(&r_surface != NULL);
    r_surface.Init();
    assert(model.GetNumberOfShellElementSets() > 0 && !model.DoShellUseMeshSurface());
    if (r_surface.GetNumberOfFacetVertices() == 3) {
      _TestXMLNodes(dynamic_cast<tledDeformableMembraneContactSurfaceT3CPU&>(r_surface));
    } else {
      std::cerr << "Incompatible contact surface.\n";
      std::abort();
    }
  }
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestXMLNodes(tledUnitTest::GetMeshPath("two_cubes_sheet.xml"));

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
