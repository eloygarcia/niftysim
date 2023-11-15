// =========================================================================
// File:       tledSurfaceTopology.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TSurface>
class tledSurfaceTopology<TSurface>::_EdgeComparator {
public:
  bool operator()(const std::pair<int, int> &e0, const std::pair<int, int> &e1) const {
    assert(e0.first < e0.second && e1.first < e1.second);

    return e0.first < e1.first || (e0.first == e1.first && e0.second < e1.second);
  }
};

template <class TSurface>
void tledSurfaceTopology<TSurface>::ComputeEdges() {
  using namespace tledVectorArithmetic;
  typedef std::map<std::pair<int, int>, int, _EdgeComparator> __EdgeMap;

  std::vector<int>::const_iterator ic_fInd;
  __EdgeMap edgeMap;
    	
  m_FacetEdges.reserve(this->GetSurface().GetNumberOfFacets());  
  for (int fInd = 0; fInd < this->GetSurface().GetNumberOfFacets(); fInd++) {
    const Facet &facet = this->GetSurface().GetFacet(fInd);

    tledArray<int, Facet::NumberOfVertices> facetEdges;

    for (int eInd = 0; eInd < TSurface::Facet::NumberOfVertices; eInd++) {
      const int v0Ind = facet.NodeIndices[eInd];
      const int v1Ind = facet.NodeIndices[(eInd+1)%TSurface::Facet::NumberOfVertices];

      std::pair<int, int> edge(std::min(v0Ind, v1Ind), std::max(v0Ind, v1Ind));
      typename __EdgeMap::iterator i_mapEntry;

      assert(edge.first != edge.second);
      assert(edge.first >= 0 && edge.first < this->GetSurface().GetNumberOfNodes());
      assert(edge.second >= 0 && edge.second < this->GetSurface().GetNumberOfNodes());
      if ((i_mapEntry = edgeMap.find(edge)) == edgeMap.end()) {
	std::pair<int, int> neigh;

	neigh.first = fInd;
	neigh.second = -1;
	facetEdges[eInd] = edgeMap[edge] = m_Edges.size();
	m_Edges.push_back(edge);	    
	m_EdgeNeighbours.push_back(neigh);
      } else {
	facetEdges[eInd] = i_mapEntry->second;      
	m_EdgeNeighbours[i_mapEntry->second].second = fInd;
      }
    } /* for facet edges */
    m_FacetEdges.push_back(facetEdges);
  } /* for facets */
  assert(m_EdgeNeighbours.size() == m_Edges.size());
  assert((int)m_FacetEdges.size() == this->GetSurface().GetNumberOfFacets());
}
