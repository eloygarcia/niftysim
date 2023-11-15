// =========================================================================
// File:       tledTrackingBVHTraverser.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBaseTraverser>
tledTrackingBVHTraverserImpl<TBaseTraverser>::tledTrackingBVHTraverserImpl(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : TBaseTraverser(r_slaveBVH, masterBVH) {
  m_TrackingStep = 0;
}

template <class TBaseTraverser>
void tledTrackingBVHTraverserImpl<TBaseTraverser>::Init(tledUnstructuredContactManager &r_manager) {
  Superclass::Init(r_manager);
  tledLogDebugStream(tledHelper::Info() << "Starting init of tracking data...");
  m_TrackingStep = 0;
  if (this->GetNumberOfTrackingIterations() < 0) {
    this->SetNumberOfTrackingIterations(this->GetSlaveMesh().GetCoordinateHistorySize()/2);
  }

  if (this->GetNumberOfTrackingIterations()%2 == 0) {
    tledLogDebugStream(tledHelper::Info() << "Even number of tracking iterations detected, decreasing to " << this->GetNumberOfTrackingIterations() - 1);
    this->SetNumberOfTrackingIterations(this->GetNumberOfTrackingIterations() - 1);
  }
}

template <class TBaseTraverser>
void tledTrackingBVHTraverserImpl<TBaseTraverser>::RunNarrowPhase() {
  Superclass::RunNarrowPhase();
  this->ExtractTrackingInformationFromContacts();
}

template <class TBaseTraverser>
void tledTrackingBVHTraverserImpl<TBaseTraverser>::RunBroadPhase() {  
  if (this->GetTrackingStep()%this->GetNumberOfTrackingIterations() == 0) {
    tledLogDebugStream(tledHelper::Info() << "Running full broad-phase...");
    this->PreFullBVHTraversalHook();
    Superclass::RunBroadPhase();
    m_TrackingStep = 0;
  } else {
    tledLogDebugStream(tledHelper::Info() << "Running tracking broad-phase..." << this->GetNumberOfTrackingIterations() - this->GetTrackingStep() - 1 << " more tracking broad phases until full contact search.");
    this->RunTrackingBroadPhase();
    m_TrackingStep += 1;
  }
}

template <class TBaseTraverser>
void tledTrackingBVHTraverserImpl<TBaseTraverser>::_ConvertNodeListsToLinearList(std::vector<int> &r_list, std::vector<std::pair<int, int> > &r_ranges, int &r_maxRangeSize, const std::vector<std::vector<int> > &nodeLists) {
  const int numNodes = nodeLists.size();

  r_ranges.clear(), r_list.clear();
  r_maxRangeSize = 0;
  r_list.reserve(numNodes);
  r_ranges.resize(numNodes);
  for (int n = 0; n < numNodes; n++) {
    r_ranges[n].first = r_list.size();
    r_list.insert(r_list.end(), nodeLists[n].begin(), nodeLists[n].end());
    r_ranges[n].second = r_list.size();
    assert(r_ranges[n].second > r_ranges[n].first && r_ranges[n].first >= 0 && r_ranges[n].second <= int(r_list.size()));
    r_maxRangeSize = std::max(r_maxRangeSize, r_ranges[n].second - r_ranges[n].first);
  }
  assert(r_ranges[0].first == 0 && r_ranges.back().second == int(r_list.size()));
}

template <class TBaseTraverser>
template <class TSurface>
void tledTrackingBVHTraverserImpl<TBaseTraverser>::ExtractNodeEdgeNeighbourHood(std::vector<int> &r_list, std::vector<std::pair<int, int> > &r_ranges, int &r_maxRangeSize, const TSurface &surface) {
  std::vector<std::vector<int> > nodeEdges;

  nodeEdges.resize(surface.GetNumberOfNodes());
  for (int e = 0; e < surface.GetNumberOfEdges(); e++) {
    nodeEdges[surface.GetEdge(e).first].push_back(e);
    nodeEdges[surface.GetEdge(e).second].push_back(e);
  }

  tledTrackingBVHTraverserImpl::_ConvertNodeListsToLinearList(r_list, r_ranges, r_maxRangeSize, nodeEdges);

#ifndef NDEBUG
  assert(int(r_ranges.size()) == surface.GetNumberOfNodes() && int(r_list.size()) > surface.GetNumberOfEdges());
  for (std::vector<std::pair<int, int> >::const_iterator ic_r = r_ranges.begin(); ic_r < r_ranges.end(); ic_r++) {
    for (int i = ic_r->first; i < ic_r->second; i++) {
      assert(r_list[i] >= 0 && r_list[i] < surface.GetNumberOfEdges());
      assert(surface.GetEdge(r_list[i]).first == int(ic_r - r_ranges.begin()) || surface.GetEdge(r_list[i]).second == int(ic_r - r_ranges.begin()));
    }
  }
#endif
}

template <class TBaseTraverser>
template <class TSurface>
void tledTrackingBVHTraverserImpl<TBaseTraverser>::ExtractNodeFacetNeighbourHood(std::vector<int> &r_list, std::vector<std::pair<int, int> > &r_ranges, int &r_maxRangeSize, const TSurface &surface) {
  const int numFacetVertices = TSurface::Facet::NumberOfVertices;

  std::vector<std::vector<int> > nodeFacets;

  nodeFacets.resize(surface.GetNumberOfNodes());
  for (int f = 0; f < surface.GetNumberOfFacets(); f++) {
    for (int v = 0; v < numFacetVertices; v++) {
      assert(surface.GetFacet(f).NodeIndices[v] < int(nodeFacets.size()) && surface.GetFacet(f).NodeIndices[v] >= 0);
      nodeFacets[surface.GetFacet(f).NodeIndices[v]].push_back(f);
    }
  }

  tledTrackingBVHTraverserImpl::_ConvertNodeListsToLinearList(r_list, r_ranges, r_maxRangeSize, nodeFacets);

#ifndef NDEBUG
  assert(int(r_ranges.size()) == surface.GetNumberOfNodes() && int(r_list.size()) > surface.GetNumberOfFacets());
  for (std::vector<std::pair<int, int> >::const_iterator ic_r = r_ranges.begin(); ic_r < r_ranges.end(); ic_r++) {
    for (int i = ic_r->first; i < ic_r->second; i++) {
      assert(int(r_list.size()) > i);
      assert(r_list[i] >= 0 && r_list[i] < surface.GetNumberOfFacets());
    }
  }
#endif
}
