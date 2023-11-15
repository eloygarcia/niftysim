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
#include "tledSelfCollisionBVH.h"
#include "tledSelfCollisionBVHCreator.h"
#include "tledSelfCollisionBVHXMLExporter.h"
#include "tledModel.h"
#include "tledSolverCPU.h"
#include "xmlParser.h"
#include "tledAABB.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>

template <class TBVH>
void _TestXMLNodes(TBVH &r_bvh) {
  tledSelfCollisionBVHXMLExporter<TBVH> exporter;
  XMLNode root, testChild;

  exporter.SetInput(r_bvh);
  exporter.Export();

  root = exporter.GetRootNode();

  tledUnitTestAssert(root.nChildNode("Margin") == 1);
  tledUnitTestAssert(root.nChildNode("LeafBVIndices") == 1);
  tledUnitTestAssert(root.nChildNode("BoundingVolumes") == 1);
  tledUnitTestAssert(root.nChildNode("NonAdjacentGeometryNodes") == 1);
  tledUnitTestAssert(root.nChildNode("GeometryClusterSubtreeRootIndices") == 1);
  tledUnitTestAssert(root.nChildNode("SelfCollisionCandidates") == 1);
  tledUnitTestAssert(root.nChildNode("BoundingVolumes") == 1);  

  testChild = root.getChildNode("BoundingVolumes");  
  tledUnitTestAssert(testChild.nChildNode("SelfCollisionBV") == r_bvh.GetNumberOfBVs());
  for (int bvInd = 0; bvInd < r_bvh.GetNumberOfBVs(); bvInd++) {
    const typename TBVH::BoundingVolume &bv = r_bvh.GetBV(bvInd);

    XMLNode bvNode = testChild.getChildNode("SelfCollisionBV", bvInd);
    std::vector<float> bounds;
    tledUnitTestAssert(bvNode.nChildNode("Bounds") == 1);
    tledUnitTestAssert(bvNode.nChildNode("ParentIndex") == 1);
    tledUnitTestAssert((bvNode.nChildNode("ChildIndices") == 1 && bvNode.nChildNode("PrimitiveIndex") == 0) || (bvNode.nChildNode("ChildIndices") == 0 && bvNode.nChildNode("PrimitiveIndex") == 1));

    bounds = GetXMLTextAsVector<float>(bvNode.getChildNode("Bounds"));
    tledUnitTestAssert(2*3 == bounds.size());
    for (int c = 0; c < 3; c++) for (int b = 0; b < 2; b++) {
	tledUnitTestAssert(std::fabs(bounds[c*2+b] - bv.Bounds[c][b]) < 1e-4f*std::fabs(bv.Bounds[c][b]));
      }
  }

  root.deleteNodeContent();
}

static void _TestXMLNodes(const std::string &inSimFile) {
  tledModel model(inSimFile.c_str());
  tledSolverCPU solver;

  solver.Init(&model);

  {
    tledUnstructuredContactManager man(model, solver, false);

    if (man.GetDeformableSurface<tledDeformableContactSurface>().GetNumberOfFacetVertices() == 3 && (model.GetNumberOfShellElementSets() == 0 || model.DoShellUseMeshSurface())) {
      _TestXMLNodes(dynamic_cast<tledSelfCollisionBVHImpl<tledDeformableContactSurfaceT3CPU, tledAABB<2> >&>(man.GetDeformableBVH()));
    } else {
      tledFatalError("Incompatible contact surface.");
    }
  }
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestXMLNodes(tledUnitTest::GetMeshPath("ball_ball.xml"));

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
