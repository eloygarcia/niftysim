// =========================================================================
// File:       tledMeshTopology.cpp
// Purpose:    Mesh topology class
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

#include "tledMeshTopology.h"

using namespace std;

template <const int t_numFacetVtcs>
struct _FacetHasher {
  size_t operator()(const typename tledMeshTopology<t_numFacetVtcs == 3? 4 : 8>::Facet &key) const {
    int cInd;
    size_t hashval;

    hashval = key[0];
    for (cInd = 1; cInd < t_numFacetVtcs; cInd++) hashval = hashval*11 + key[cInd];

    return hashval;
  }
};

template <const int t_numFacetNodes, const int t_numElementNodes>
static void _ExtractFacets(std::vector<typename tledMeshTopology<t_numElementNodes>::Facet> &r_facets, std::vector<typename tledMeshTopology<t_numElementNodes>::ElementFacetList> &r_elementFacetMap, std::vector<std::vector<int> > &r_nodeFacetMap,
			   const int localFacetDefs[][t_numFacetNodes], const int elDefs[], const int numEls, const int numNodes) {
  typedef typename tledMeshTopology<t_numElementNodes>::Facet __Facet;
  typedef tledMap<__Facet, int, _FacetHasher<t_numFacetNodes> > __FacetMap;

  const int numElFacets = t_numElementNodes == 4? 4 : 6;

  __FacetMap facetMap;
  std::vector<int> numValidElFacets;

  numValidElFacets.insert(numValidElFacets.end(), numEls, 0);
  r_nodeFacetMap.resize(numNodes);
  r_elementFacetMap.resize(numEls);
  r_facets.reserve(numEls/2);

  for (int eInd = 0; eInd < numEls; eInd++) {
    const int *elnodes = elDefs + t_numElementNodes*eInd;

    for (int elFacetInd = 0; elFacetInd < numElFacets; elFacetInd++) {	
      typename __FacetMap::const_iterator ic_facet;
      __Facet newFacet(elnodes[localFacetDefs[elFacetInd][0]], elnodes[localFacetDefs[elFacetInd][1]], elnodes[localFacetDefs[elFacetInd][2]], t_numFacetNodes == 4? elnodes[localFacetDefs[elFacetInd][3]] : -1);

      if ((ic_facet = facetMap.find(newFacet)) == facetMap.end()) {
	const int facetInd = r_facets.size();

	int vInd;

	r_facets.push_back(newFacet);
	facetMap.insert(typename __FacetMap::value_type(newFacet, facetInd));
	r_elementFacetMap[eInd][numValidElFacets[eInd]] = facetInd;
	for (vInd = 0; vInd < t_numFacetNodes; vInd++) r_nodeFacetMap[newFacet[vInd]].push_back(facetInd);

	assert((*facetMap.find(newFacet)).second == facetInd);
      } else r_elementFacetMap[eInd][numValidElFacets[eInd]] = ic_facet->second;
	
      numValidElFacets[eInd] += 1;
    }
  }  
}

template <const int t_numFacetNodes, const int t_numElFacets>
static void _ExtractFacetNeighbours(std::vector<std::pair<int, int> > &r_facetEls, const std::vector<tledArray<int, t_numElFacets> > &elementFacetMap, const int numFacets) {
  static const std::pair<int,int> undefined(-1, -1);

  const int numEls = elementFacetMap.size();

  r_facetEls.clear();
  r_facetEls.insert(r_facetEls.end(), numFacets, undefined);
  for (int elInd = 0; elInd < numEls; elInd++) {
    const tledArray<int, t_numElFacets> &elemfacets = elementFacetMap[elInd];
    
    for (int elFacetInd = 0; elFacetInd < t_numElFacets; elFacetInd++) {		
      if (r_facetEls[elemfacets[elFacetInd]].first == -1) r_facetEls[elemfacets[elFacetInd]].first = elInd;
      else {
	assert(r_facetEls[elemfacets[elFacetInd]].second == -1);
	r_facetEls[elemfacets[elFacetInd]].second = elInd;
      }
    }
  }
}

template <>
void tledMeshTopology<8>::ComputeFacets() {
  static const int hexFacets[][4] = {{0, 1, 2, 3}, {4, 5, 6, 7},
				     {0, 1, 5, 4}, {1, 2, 6, 5},
				     {2, 3, 7, 6}, {3, 0, 4, 7}};

  _ExtractFacets<4, 8>(m_Facets, m_ElFacets, m_NodeFacets, hexFacets, GetMesh().GetAllElNodeInds(), GetMesh().GetNumEls(), GetMesh().GetNumNodes());
  _ExtractFacetNeighbours<4, 6>(m_FacetEls, m_ElFacets, GetNumFacets());
} 

template <>
void tledMeshTopology<4>::ComputeFacets() {
  static const int tetraFacets[][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};				       

  _ExtractFacets<3, 4>(m_Facets, m_ElFacets, m_NodeFacets, tetraFacets, GetMesh().GetAllElNodeInds(), GetMesh().GetNumEls(), GetMesh().GetNumNodes());
  _ExtractFacetNeighbours<3, 4>(m_FacetEls, m_ElFacets, GetNumFacets());
} 
