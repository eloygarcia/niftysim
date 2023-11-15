// =========================================================================
// File:       testMeshTopology.cpp
// Purpose:    Mesh topology module unit test
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
#include "tledModel.h"
#include "tledDeformableContactSurfaceCPU.h"

#include <cstdlib>
#include <string>

/*
 * Tests that all node facets are indeed adjacent to the node and that there are no adjacent facets being missed.
 */
static void _TestEdges(const tledMesh &mesh) {
  tledDeformableContactSurfaceT3CPU surface(mesh);
  std::vector<std::pair<int, int> >::const_iterator ic_edge;
  std::vector<std::pair<int, int> > neighbours;
  std::vector<std::pair<int, int> >::const_iterator ic_neighbours;
  
  surface.Init();

  neighbours.reserve(surface.GetAllEdges().size());
  for (ic_edge = surface.GetAllEdges().begin(); ic_edge < surface.GetAllEdges().end(); ic_edge++) {
    std::vector<tledDeformableContactSurfaceImpl<3>::Facet>::const_iterator ic_facet;

    tledUnitTestAssert(ic_edge->first >= 0 && ic_edge->first < surface.GetNumberOfNodes());
    tledUnitTestAssert(ic_edge->second >= 0 && ic_edge->second < surface.GetNumberOfNodes());

    neighbours.push_back(std::pair<int, int>(-1, -1));
    for (ic_facet = surface.GetAllFacets().begin(); ic_facet < surface.GetAllFacets().end(); ic_facet++) {
      if ((ic_edge->first == ic_facet->NodeIndices[0] || ic_edge->first == ic_facet->NodeIndices[1] || ic_edge->first == ic_facet->NodeIndices[2])
	  && (ic_edge->second == ic_facet->NodeIndices[0] || ic_edge->second == ic_facet->NodeIndices[1] || ic_edge->second == ic_facet->NodeIndices[2])) {
	if (neighbours.back().first == -1) neighbours.back().first = ic_facet - surface.GetAllFacets().begin();
	else if (neighbours.back().second == -1) neighbours.back().second = ic_facet - surface.GetAllFacets().begin();
	else {
	  tledUnitTestAssert(neighbours.back().second == -1 || neighbours.back().first == -1);
	}
      }
    }
  }

  for (ic_neighbours = neighbours.begin(); ic_neighbours < neighbours.end(); ic_neighbours++) {
    tledUnitTestAssert(ic_neighbours->first >= 0 && ic_neighbours->first < surface.GetNumberOfFacets());
    tledUnitTestAssert(ic_neighbours->first >= 0 && ic_neighbours->first < surface.GetNumberOfFacets());

    {
      const int edgeInd = ic_neighbours - neighbours.begin();
      const tledDeformableContactSurfaceImpl<3>::Facet &facet0 = surface.GetFacet(ic_neighbours->first);
      const tledDeformableContactSurfaceImpl<3>::Facet &facet1 = surface.GetFacet(ic_neighbours->second);

      tledUnitTestAssert(facet0.EdgeIndices[0] == edgeInd || facet0.EdgeIndices[1] == edgeInd || facet0.EdgeIndices[2] == edgeInd);
      tledUnitTestAssert(facet1.EdgeIndices[0] == edgeInd || facet1.EdgeIndices[1] == edgeInd || facet1.EdgeIndices[2] == edgeInd);
    }
  }
} /* _TestEdges */

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestEdges(tledUnitTest::LoadMSHMesh(tledUnitTest::GetMeshPath("collision_bar.msh"), "T4"));
  _TestEdges(tledUnitTest::LoadMSHMesh(tledUnitTest::GetMeshPath("organic_shape.msh"), "T4"));
  _TestEdges(tledUnitTest::LoadMSHMesh(tledUnitTest::GetMeshPath("sphere.msh"), "T4"));

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
