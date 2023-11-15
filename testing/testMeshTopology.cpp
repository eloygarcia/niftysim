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
#include "tledMeshTopology.h"

#include <cstdlib>

using namespace std;

/*
 * Tests that all node facets are indeed adjacent to the node and that there are no adjacent facets being missed.
 */
template <const int t_numElNodes>
static void _TestNodeFacets(const tledMeshTopology<t_numElNodes> &topo, const tledMesh &mesh) {
  const int numFacetNodes = t_numElNodes == 4? 3 : 4;

  for (int nodeInd = 0; nodeInd < mesh.GetNumNodes(); nodeInd++) {
    const std::vector<int> &nodeFacetInds = topo.GetNodeFacets(nodeInd);

    std::vector<int>::const_iterator ic_facetInd;
    int meshFacetInd;

    for (ic_facetInd = nodeFacetInds.begin(); ic_facetInd < nodeFacetInds.end(); ic_facetInd++) {
      const typename tledMeshTopology<t_numElNodes>::Facet &facet = topo.GetFacet(*ic_facetInd);

      tledUnitTestAssert(facet[0] == nodeInd || facet[1] == nodeInd || facet[2] == nodeInd || (numFacetNodes == 4 && facet[3] == nodeInd));
    }

    for (meshFacetInd = 0; meshFacetInd < topo.GetNumFacets(); meshFacetInd++) {
      const typename tledMeshTopology<t_numElNodes>::Facet &facet = topo.GetFacet(meshFacetInd);

      if (facet[0] == nodeInd || facet[1] == nodeInd || facet[2] == nodeInd || (numFacetNodes == 4 && facet[3] == nodeInd)) {
	for (ic_facetInd = nodeFacetInds.begin(); ic_facetInd < nodeFacetInds.end() && *ic_facetInd != meshFacetInd; ic_facetInd++);
	tledUnitTestAssert(ic_facetInd < nodeFacetInds.end());
      }
    }
  }
} /* _TestNodeFacets */

/**
 * \brief Runs a number of tests on the given geometry
 */
template <const int t_numElNodes>
static void _TestTopology(const string &meshPath) {
  tledMesh mesh = tledUnitTest::LoadMSHMesh(meshPath, t_numElNodes == 4? "T4" : "H8");

  {
    tledMeshTopology<t_numElNodes> topo(mesh);

    topo.ComputeNodeAdjacency();
    topo.ComputeFacets();
    
    _TestNodeFacets<t_numElNodes>(topo, mesh);
  } 
} /* _TestTopology */

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestTopology<8>(tledUnitTest::GetMeshPath("box8.msh"));

  _TestTopology<4>(tledUnitTest::GetMeshPath("collision_bar.msh"));
  _TestTopology<4>(tledUnitTest::GetMeshPath("organic_shape.msh"));
  _TestTopology<4>(tledUnitTest::GetMeshPath("sphere.msh"));

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
